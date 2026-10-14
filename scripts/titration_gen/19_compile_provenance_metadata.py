#!/usr/bin/env python3
"""
Compile a single, citable provenance record for the bacterial genome panel,
human reads, and viral titration dataset used across this project's HRR
benchmarking.

Writes:
  metadata/bacterial_genomes.tsv      - one row per genome, + source/publication
  metadata/human_reads.tsv            - one row per Forbes-27 ENA run, + source/publication
  metadata/human_reads_african25.tsv  - one row per continental-African ENA run (5/site x 5 sites)
  metadata/viral_reads.tsv           - one row per *actually used* viral ENA run
                                        (the exact 27 accessions in
                                        _ena_utils.BENCHMARK_VIRUS_TYPES x
                                        BENCHMARK_ITERATION - not the full
                                        351-accession catalog Forbes et al.
                                        deposited), joined from their
                                        titration_sheet.csv (Sample/titre/
                                        iteration/study/sample accession) and
                                        err_to_titre.tsv (adds the run
                                        accession) on ena_sample_accession
"""
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _ena_utils import BENCHMARK_ITERATION, BENCHMARK_VIRUS_TYPES  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
METADATA_DIR = REPO_ROOT / "metadata"

GENOME_MANIFEST_PATH = REPO_ROOT / "results" / "genome_manifest.tsv"
HUMAN_ACCESSIONS_PATH = REPO_ROOT / "results" / "human_accessions.tsv"
HUMAN_ACCESSIONS_AFRICAN25_PATH = REPO_ROOT / "results" / "human_accessions_african25.tsv"

# Forbes et al. (2025) sample sheets for ENA project ERP169281 (not included in
# this repository; metadata/viral_reads.tsv is the output already provided).
VIRAL_TITRATION_SHEET_PATH = Path(os.environ.get("VIRAL_TITRATION_SHEET", REPO_ROOT / "data" / "titration_sheet.csv"))
VIRAL_ERR_TO_TITRE_PATH = Path(os.environ.get("VIRAL_ERR_TO_TITRE", REPO_ROOT / "data" / "err_to_titre.tsv"))

BACTERIAL_SOURCE_DB = "NCBI RefSeq (https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/)"
BACTERIAL_SELECTION_PUBLICATION = (
    "Maghini DG, Oduaran OH, Olubayo LAI, et al. Expanding the human gut microbiome atlas of "
    "Africa. Nature 2025;638:718-728. https://doi.org/10.1038/s41586-024-08485-8 "
    "(Supplementary Data 3 species-abundance ranking)"
)

HUMAN_SOURCE_DB = "ENA / 1000 Genomes Project 30x WGS (study PRJEB31736)"
HUMAN_DATA_PUBLICATION = (
    "1000 Genomes Project Consortium. A global reference for human genetic variation. "
    "Nature 2015;526:68-74. https://doi.org/10.1038/nature15393"
)
HUMAN_DESIGN_PUBLICATION = (
    "Forbes M, Ng DYK, Boggan RM, et al. Benchmarking of human read removal strategies for "
    "viral and microbial metagenomics. Cell Reports Methods 2025;5:101218. "
    "https://doi.org/10.1016/j.crmeth.2025.101218 (mixed-ancestry population/accession design)"
)

VIRAL_SOURCE_DB = "ENA project ERP169281 (Forbes et al. deposited clinical virome + human titration dataset)"
VIRAL_PUBLICATION = (
    "Forbes M, Ng DYK, Boggan RM, et al. Benchmarking of human read removal strategies for "
    "viral and microbial metagenomics. Cell Reports Methods 2025;5:101218. "
    "https://doi.org/10.1016/j.crmeth.2025.101218"
)

SUPERPOPULATION = {
    "CEU": "EUR", "TSI": "EUR", "GBR": "EUR", "FIN": "EUR", "IBS": "EUR",
    "YRI": "AFR", "LWK": "AFR", "ASW": "AFR", "ACB": "AFR", "MSL": "AFR", "ESN": "AFR", "GWD": "AFR",
    "CHB": "EAS", "JPT": "EAS", "CHS": "EAS", "KHV": "EAS", "CDX": "EAS",
    "GIH": "SAS", "PJL": "SAS", "BEB": "SAS", "STU": "SAS", "ITU": "SAS",
    "MXL": "AMR", "PUR": "AMR", "CLM": "AMR", "PEL": "AMR",
}


def compile_bacterial():
    with open(GENOME_MANIFEST_PATH, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    out_rows = []
    for r in rows:
        out_rows.append(
            {
                "binomial": r["binomial"],
                "refseq_accession": r["genome_ID"],
                "taxid": r["taxid"],
                "species_taxid": r["species_taxid"],
                "awi_gen_abundance_rank": r["rank"],
                "awi_gen_abundance_pct": r["abundance_pct"],
                "assembly_level": r["assembly_level"],
                "source_database": BACTERIAL_SOURCE_DB,
                "selection_rationale_publication": BACTERIAL_SELECTION_PUBLICATION,
            }
        )
    return out_rows


def compile_human(accessions_path):
    with open(accessions_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    out_rows = []
    for r in rows:
        out_rows.append(
            {
                "population_code": r["population_code"],
                "population_name": r["population_name"],
                "superpopulation": SUPERPOPULATION.get(r["population_code"], ""),
                "ena_run_accession": r["run_accession"],
                "ena_study_accession": r["ena_study_accession"],
                "source_database": HUMAN_SOURCE_DB,
                "data_publication": HUMAN_DATA_PUBLICATION,
                "design_rationale_publication": HUMAN_DESIGN_PUBLICATION,
            }
        )
    return out_rows


def compile_viral():
    with open(VIRAL_TITRATION_SHEET_PATH, newline="") as f:
        sample_rows = list(csv.DictReader(f))
    with open(VIRAL_ERR_TO_TITRE_PATH, newline="") as f:
        run_rows = list(csv.DictReader(f, delimiter="\t"))

    run_by_sample_acc = {r["ena_sample_accession"]: r["ena_run_accession"] for r in run_rows}

    out_rows = []
    missing = 0
    for r in sample_rows:
        if r["Sample"] not in BENCHMARK_VIRUS_TYPES or int(r["iteration"]) != BENCHMARK_ITERATION:
            continue  # only the accessions actually used (20_extract_viral_pool.py), not the full catalog
        sample_acc = r["ena_sample_accession"]
        run_acc = run_by_sample_acc.get(sample_acc)
        if run_acc is None:
            missing += 1
            continue
        out_rows.append(
            {
                "virus_sample_type": r["Sample"],
                "human_titre_decile": r["human_titre"],
                "iteration": r["iteration"],
                "ena_study_accession": r["ena_study_accession"],
                "ena_sample_accession": sample_acc,
                "ena_run_accession": run_acc,
                "source_database": VIRAL_SOURCE_DB,
                "publication": VIRAL_PUBLICATION,
            }
        )
    if missing:
        print(f"WARNING: {missing} titration_sheet.csv rows had no matching run accession in err_to_titre.tsv")
    return out_rows


def write_tsv(path, rows):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {path}")


def main():
    METADATA_DIR.mkdir(parents=True, exist_ok=True)

    write_tsv(METADATA_DIR / "bacterial_genomes.tsv", compile_bacterial())
    write_tsv(METADATA_DIR / "human_reads.tsv", compile_human(HUMAN_ACCESSIONS_PATH))
    write_tsv(METADATA_DIR / "human_reads_african25.tsv", compile_human(HUMAN_ACCESSIONS_AFRICAN25_PATH))
    write_tsv(METADATA_DIR / "viral_reads.tsv", compile_viral())


if __name__ == "__main__":
    main()
