#!/usr/bin/env python3
"""
Resolve 5 individuals per continental-African site (25 total: YRI, LWK,
GWD, MSL, ESN x 5 each) instead of 1, so the African-ancestry human
background isn't built from a single, potentially idiosyncratic genome per
site - averaging across 5 independent individuals per population reduces
that single-genome noise.

Two live sources are joined (neither hardcoded, both re-fetched fresh):
  1. The 1000 Genomes/IGSR pedigree file (individual_id -> population) -
     the canonical sample/population mapping, NOT the same as this repo's
     earlier per-row "ena_study_accession" column (which, on inspection,
     didn't actually resolve to anything population-specific - confirmed
     by directly querying ENA for one of our own run accessions, whose
     real study_accession is just the PRJEB31736 umbrella project, not the
     ERP1143xx value originally in data/forbes_human_ancestry_accessions.tsv).
  2. PRJEB31736's own run list - sample_title is consistently
     "Coriell GM<digits>" (-> individual NA<digits>) or
     "Coriell HG<digits>" (-> individual HG<digits>) for all 37,090 runs,
     confirmed by checking every row, not assumed.

Selection: each population's already-used individual (from
data/forbes_human_ancestry_accessions.tsv) is kept for continuity with
already-downloaded data, plus 4 more chosen via a seeded random sample
(same repo-wide seed) from the remaining matched individuals. One run
accession per individual (lowest accession, for determinism, in case an
individual has multiple sequencing-lane runs).

Output: data/continental_african_5x5_accessions.tsv - same column schema
as data/forbes_human_ancestry_accessions.tsv (population_code,
population_name, ena_run_accession, ena_study_accession), so
15_download_human_reads.py consumes it unchanged.
"""
import csv
import random
import re
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
EXISTING_ACCESSIONS_PATH = REPO_ROOT / "data" / "forbes_human_ancestry_accessions.tsv"
OUTPUT_PATH = REPO_ROOT / "data" / "continental_african_5x5_accessions.tsv"

PED_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped"
PRJEB31736_PROJECT = "PRJEB31736"

TARGET_POPULATIONS = {"YRI", "LWK", "GWD", "MSL", "ESN"}
INDIVIDUALS_PER_POPULATION = 5
SELECTION_SEED = 632741178  # same repo-wide seed used elsewhere

POPULATION_NAMES = {
    "YRI": "Yoruba in Ibadan, Nigeria",
    "LWK": "Luhya in Webuye, Kenya",
    "GWD": "Gambian in Western Division, The Gambia",
    "MSL": "Mende in Sierra Leone",
    "ESN": "Esan in Nigeria",
}

CORIELL_TITLE_RE = re.compile(r"Coriell (GM|HG)(\d+)")


def load_existing_individual_per_population():
    """The 1 already-used run accession per population - kept for
    continuity. Returns {population_code: ena_run_accession}."""
    with open(EXISTING_ACCESSIONS_PATH, newline="") as f:
        rows = {r["population_code"]: r["ena_run_accession"] for r in csv.DictReader(f, delimiter="\t")}
    return {pop: rows[pop] for pop in TARGET_POPULATIONS}


def fetch_ped_population_map():
    with urllib.request.urlopen(PED_URL, timeout=60) as resp:
        text = resp.read().decode("utf-8")
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    return {row["Individual ID"]: row["Population"] for row in reader if row["Population"] in TARGET_POPULATIONS}


def fetch_prjeb31736_runs():
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={PRJEB31736_PROJECT}&result=read_run&fields=run_accession,sample_accession,sample_title&format=tsv"
    )
    with urllib.request.urlopen(url, timeout=120) as resp:
        text = resp.read().decode("utf-8")
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def individual_id_from_title(sample_title):
    m = CORIELL_TITLE_RE.match(sample_title)
    if not m:
        return None
    prefix, digits = m.groups()
    return f"NA{digits}" if prefix == "GM" else f"HG{digits}"


def main():
    existing_by_pop = load_existing_individual_per_population()

    print(f"Fetching 1000 Genomes pedigree from {PED_URL} ...")
    ped_pop = fetch_ped_population_map()
    print(f"{len(ped_pop)} individuals across {TARGET_POPULATIONS} in the pedigree file.")

    print(f"Fetching {PRJEB31736_PROJECT} run list ...")
    runs = fetch_prjeb31736_runs()
    print(f"{len(runs)} runs in {PRJEB31736_PROJECT}.")

    runs_by_individual = {}
    for row in runs:
        indiv = individual_id_from_title(row["sample_title"])
        if indiv is None:
            continue
        runs_by_individual.setdefault(indiv, []).append(row["run_accession"])

    matched_by_pop = {pop: [] for pop in TARGET_POPULATIONS}
    for indiv, pop in ped_pop.items():
        if indiv in runs_by_individual:
            matched_by_pop[pop].append(indiv)

    out_rows = []
    for pop in sorted(TARGET_POPULATIONS):
        candidates = sorted(matched_by_pop[pop])
        existing_run = existing_by_pop[pop]
        rng = random.Random(SELECTION_SEED + hash(pop) % 1000)

        # exclude whichever individual the existing run already belongs to,
        # so the random pick can't accidentally duplicate it
        existing_individual = next(
            (indiv for indiv in candidates if existing_run in runs_by_individual[indiv]), None
        )
        pool = [i for i in candidates if i != existing_individual]
        chosen_new = rng.sample(pool, INDIVIDUALS_PER_POPULATION - 1)
        chosen_individuals = ([existing_individual] if existing_individual else []) + chosen_new
        # top up if the existing run's individual wasn't found in the matched set
        while len(chosen_individuals) < INDIVIDUALS_PER_POPULATION:
            extra = rng.choice([i for i in pool if i not in chosen_individuals])
            chosen_individuals.append(extra)

        print(f"{pop}: {len(candidates)} candidates, chose {chosen_individuals}")
        for indiv in chosen_individuals:
            run_acc = existing_run if indiv == existing_individual else sorted(runs_by_individual[indiv])[0]
            out_rows.append(
                {
                    "population_code": pop,
                    "population_name": POPULATION_NAMES[pop],
                    "individual_id": indiv,
                    "ena_run_accession": run_acc,
                    "ena_study_accession": PRJEB31736_PROJECT,
                }
            )

    with open(OUTPUT_PATH, "w", newline="") as f:
        fieldnames = ["population_code", "population_name", "ena_run_accession", "ena_study_accession", "individual_id"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in out_rows:
            writer.writerow({k: r[k] for k in fieldnames})

    print(f"Wrote {len(out_rows)} rows ({INDIVIDUALS_PER_POPULATION} individuals x {len(TARGET_POPULATIONS)} populations) to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
