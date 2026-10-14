#!/usr/bin/env python3
"""
Resolve and download a bacterial RefSeq genome panel grounded in the real
measured species-abundance ranking of the AWI-Gen 2 cohort (Maghini,
Oduaran, Olubayo et al. 2025, Nature 638:718-728, Supplementary Data 3),
instead of a hand-curated "known gut taxa" list.

data/top60_candidates.csv (60 GTDB-named mOTU entries, ranked by aggregated
mean relative abundance across the 6 AWI-Gen 2 study sites) is parsed into
deduplicated species candidates (GTDB splits some species into per-clade
suffixed names, e.g. Faecalibacterium prausnitzii_D/_C/_G/_E, which are the
SAME named species and must collapse to one), then every candidate -
including GTDB placeholder/uncultured taxa - is attempted against NCBI's
live assembly_summary.txt. Every attempt (success or failure, with a
specific reason) is logged to results/species_resolution_log.tsv; only
successes go into results/genome_manifest.tsv.
"""
import csv
import gzip
import io
import re
import shutil
import sys
import urllib.request
from pathlib import Path

ASSEMBLY_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
CANDIDATES_CSV_PATH = REPO_ROOT / "data" / "top60_candidates.csv"
GENOMES_DIR = REPO_ROOT / "genomes"
MANIFEST_PATH = REPO_ROOT / "results" / "genome_manifest.tsv"
RESOLUTION_LOG_PATH = REPO_ROOT / "results" / "species_resolution_log.tsv"

REFERENCE_CATEGORIES = {"reference genome", "representative genome"}
ASSEMBLY_LEVEL_RANK = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}

MOTU_TAG_RE = re.compile(r"\s*\[[^\]]*\]\s*$")
GTDB_SUFFIX_RE = re.compile(r"_[A-Z]+\d*$")
PLACEHOLDER_GENUS_RE = re.compile(r"^[A-Z0-9-]+$")
PLACEHOLDER_SPECIES_RE = re.compile(r"^sp\d+$")
NO_SPECIES_LEVEL_MARKERS = ("incertae sedis", "not_annotated")

# Documented recent NCBI genus reclassifications relevant to this candidate
# list, not derivable from the CSV's own columns (GTDB and NCBI taxonomy
# diverge here): Bacteroides dorei/vulgatus -> Phocaeicola (2022); Prevotella
# copri -> Segatella (2023). Confirmed against current assembly_summary.txt
# organism_name values, not assumed.
GENUS_RECLASSIFICATIONS = {
    "Bacteroides": ["Phocaeicola"],
    "Prevotella": ["Segatella"],
}


def expand_genus_aliases(genus, species):
    pairs = [(genus, species)] + [(alt, species) for alt in GENUS_RECLASSIFICATIONS.get(genus, [])]
    return [f"{g} {s}" for g, s in pairs]


def strip_motu_tag(name):
    return MOTU_TAG_RE.sub("", name).strip()


def is_no_species_level(name):
    lname = name.lower()
    return any(marker in lname for marker in NO_SPECIES_LEVEL_MARKERS)


def normalize_binomial(name):
    """('Genus_X species_Y', ...) -> ('Genus', 'species') with GTDB clade
    suffixes stripped, or None if not a clean 2-token binomial."""
    tokens = name.split(" ", 1)
    if len(tokens) != 2:
        return None
    genus, species = tokens
    return GTDB_SUFFIX_RE.sub("", genus), GTDB_SUFFIX_RE.sub("", species)


def is_placeholder_genus(genus_norm):
    return bool(PLACEHOLDER_GENUS_RE.match(genus_norm.strip("[]")))


def is_placeholder_species(species_norm):
    return bool(PLACEHOLDER_SPECIES_RE.match(species_norm))


def build_candidates(csv_path):
    """Parse top60_candidates.csv (already rank-sorted ascending) into
    deduplicated candidates, keyed by normalized binomial."""
    with open(csv_path, newline="") as f:
        rows = list(csv.DictReader(f))

    candidates = {}
    order = []

    for row in rows:
        rank = int(row["rank"])
        gtdb_name = row["gtdb_name"]
        abundance_pct = float(row["abundance_pct"])
        prevalence = float(row["prevalence"])
        original_example = row["original_species_examples"]

        bare = strip_motu_tag(gtdb_name)

        if is_no_species_level(bare):
            key = f"__no_species_level_rank{rank}__"
            candidates[key] = {
                "ranks": [rank], "gtdb_names": [gtdb_name], "binomial": None,
                "abundance_pct": abundance_pct, "prevalence": prevalence,
                "aliases": [], "no_species_level": True,
                "placeholder_genus": False, "placeholder_species": False,
            }
            order.append(key)
            continue

        norm = normalize_binomial(bare)
        genus_norm, species_norm = norm if norm else (bare, "")
        key = f"{genus_norm} {species_norm}".strip()

        if key in candidates:
            candidates[key]["ranks"].append(rank)
            candidates[key]["gtdb_names"].append(gtdb_name)
            continue  # CSV is rank-ascending, so first-seen is already best-ranked

        species_alternatives = species_norm.split("/") if "/" in species_norm else [species_norm]
        aliases = []
        for sp in species_alternatives:
            aliases.extend(expand_genus_aliases(genus_norm, sp))
        if bare != key:
            aliases.append(bare)  # literal un-stripped form, cheap to also try

        if original_example and original_example != gtdb_name:
            example_bare = strip_motu_tag(original_example)
            if not is_no_species_level(example_bare):
                example_norm = normalize_binomial(example_bare)
                if example_norm:
                    ex_genus, ex_species = example_norm
                    # only a valid synonym alias if it's a rename (same
                    # species epithet, different genus) - never substitute a
                    # genuinely different species
                    if ex_species == species_norm and ex_genus != genus_norm:
                        aliases.append(f"{ex_genus} {ex_species}")
                        if example_bare != f"{ex_genus} {ex_species}":
                            aliases.append(example_bare)

        candidates[key] = {
            "ranks": [rank], "gtdb_names": [gtdb_name], "binomial": key,
            "abundance_pct": abundance_pct, "prevalence": prevalence,
            "aliases": aliases, "no_species_level": False,
            "placeholder_genus": is_placeholder_genus(genus_norm),
            "placeholder_species": is_placeholder_species(species_norm),
        }
        order.append(key)

    return [candidates[k] for k in order]


def stream_assembly_summary_rows():
    req = urllib.request.Request(ASSEMBLY_SUMMARY_URL, headers={"User-Agent": "curl/8"})
    with urllib.request.urlopen(req) as resp:
        text_stream = io.TextIOWrapper(resp, encoding="utf-8")
        header = None
        for line in text_stream:
            line = line.rstrip("\n")
            if line.startswith("#assembly_accession"):
                header = line[1:].split("\t")
                continue
            if line.startswith("#") or not line:
                continue
            if header is None:
                continue
            yield dict(zip(header, line.split("\t")))


def organism_matches(organism_name, aliases):
    for candidate in aliases:
        if organism_name == candidate or organism_name.startswith(candidate + " "):
            return True
    return False


def select_genomes(rows, candidates):
    attemptable = [(i, c) for i, c in enumerate(candidates) if not c["no_species_level"]]
    buckets = {i: [] for i, _ in attemptable}

    for row in rows:
        if row.get("assembly_level") not in ASSEMBLY_LEVEL_RANK:
            continue
        if row.get("version_status") != "latest":
            continue
        organism = row.get("organism_name", "")
        for i, c in attemptable:
            if organism_matches(organism, c["aliases"]):
                buckets[i].append(row)
                break

    manifest_entries = []
    log_rows = []
    for i, c in enumerate(candidates):
        common = {
            "ranks": ",".join(str(r) for r in sorted(c["ranks"])),
            "gtdb_names": "; ".join(c["gtdb_names"]),
            "binomial": c["binomial"] or "",
            "abundance_pct": c["abundance_pct"],
            "prevalence": c["prevalence"],
            "aliases_tried": "; ".join(c["aliases"]),
        }

        if c["no_species_level"]:
            log_rows.append({**common, "outcome": "excluded", "accession": "",
                              "matched_organism_name": "",
                              "exclusion_reason": "no species-level classification (GTDB placeholder for unassigned reads)"})
            continue

        bucket = buckets[i]
        if not bucket:
            if c["placeholder_genus"] or c["placeholder_species"]:
                reason = "no RefSeq assembly found for GTDB placeholder taxon, likely MAG-only representation"
            else:
                reason = f"no RefSeq assembly found under any tried alias ({common['aliases_tried']})"
            log_rows.append({**common, "outcome": "excluded", "accession": "",
                              "matched_organism_name": "", "exclusion_reason": reason})
            continue

        bucket.sort(key=lambda r: r.get("seq_rel_date", ""), reverse=True)
        bucket.sort(
            key=lambda r: (
                ASSEMBLY_LEVEL_RANK[r["assembly_level"]],
                0 if r.get("refseq_category") in REFERENCE_CATEGORIES else 1,
            )
        )
        best = bucket[0]
        if best["assembly_level"] != "Complete Genome":
            print(f"NOTE: {c['binomial']} has no Complete Genome assembly, using best available level: {best['assembly_level']}", file=sys.stderr)
        manifest_entries.append((c, best))
        log_rows.append({**common, "outcome": "included", "accession": best["assembly_accession"],
                          "matched_organism_name": best["organism_name"], "exclusion_reason": ""})

    return manifest_entries, log_rows


def download_genome(accession, ftp_path, dest_fna):
    asm_name = ftp_path.rstrip("/").split("/")[-1]
    url = f"{ftp_path}/{asm_name}_genomic.fna.gz"
    gz_path = dest_fna.with_suffix(dest_fna.suffix + ".gz")
    req = urllib.request.Request(url, headers={"User-Agent": "curl/8"})
    with urllib.request.urlopen(req) as resp, open(gz_path, "wb") as out:
        shutil.copyfileobj(resp, out)
    with gzip.open(gz_path, "rb") as gz_in, open(dest_fna, "wb") as fna_out:
        shutil.copyfileobj(gz_in, fna_out)
    gz_path.unlink()


def main():
    GENOMES_DIR.mkdir(parents=True, exist_ok=True)
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)

    candidates = build_candidates(CANDIDATES_CSV_PATH)
    print(f"Parsed {len(candidates)} deduplicated candidates from {CANDIDATES_CSV_PATH.name}.")

    print("Streaming NCBI RefSeq bacterial assembly_summary.txt ...")
    rows = list(stream_assembly_summary_rows())
    print(f"Parsed {len(rows)} assembly_summary rows.")

    manifest_entries, log_rows = select_genomes(rows, candidates)
    print(f"Resolved {len(manifest_entries)} of {len(candidates)} candidates to a RefSeq assembly.")

    manifest_rows = []
    for c, row in manifest_entries:
        accession = row["assembly_accession"]
        dest_fna = GENOMES_DIR / f"{accession}.fna"
        if dest_fna.exists():
            print(f"Already have {accession} ({c['binomial']}), skipping download.")
        else:
            print(f"Downloading {accession} ({c['binomial']}) ...")
            download_genome(accession, row["ftp_path"], dest_fna)
        manifest_rows.append(
            {
                "genome_ID": accession,
                "binomial": c["binomial"],
                "rank": ",".join(str(r) for r in sorted(c["ranks"])),
                "abundance_pct": c["abundance_pct"],
                "prevalence": c["prevalence"],
                "taxid": row["taxid"],
                "species_taxid": row["species_taxid"],
                "refseq_category": row["refseq_category"],
                "assembly_level": row["assembly_level"],
                "seq_rel_date": row["seq_rel_date"],
                "ftp_path": row["ftp_path"],
                "local_path": str(dest_fna.relative_to(REPO_ROOT)),
            }
        )

    with open(MANIFEST_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(manifest_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Wrote manifest for {len(manifest_rows)} genomes to {MANIFEST_PATH}")

    with open(RESOLUTION_LOG_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(log_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(log_rows)
    n_excluded = sum(1 for r in log_rows if r["outcome"] == "excluded")
    print(f"Wrote resolution log for {len(log_rows)} candidates ({n_excluded} excluded) to {RESOLUTION_LOG_PATH}")


if __name__ == "__main__":
    main()
