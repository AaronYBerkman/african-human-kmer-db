#!/usr/bin/env python3
"""
Build CAMISIM's two required input files from genome_manifest.tsv:

  genome_locations.tsv:  genome_ID \t absolute_path_to_fasta
  metadata.tsv:          genome_ID \t OTU \t NCBI_ID \t novelty_category

OTU is set to each genome's species_taxid (so the 4 deliberately-duplicated
species share an OTU, letting CAMISIM's min/max_strains_per_otu actually do
strain-level sampling) and NCBI_ID is the genome's own (more specific) taxid.
All genomes are named, classified RefSeq species, so novelty_category is
"known_strain" throughout.
"""
import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
MANIFEST_PATH = REPO_ROOT / "results" / "genome_manifest.tsv"
CAMISIM_DIR = REPO_ROOT / "camisim_inputs"
GENOME_LOCATIONS_PATH = CAMISIM_DIR / "genome_locations.tsv"
METADATA_PATH = CAMISIM_DIR / "metadata.tsv"


def main():
    with open(MANIFEST_PATH, newline="") as f:
        genomes = list(csv.DictReader(f, delimiter="\t"))

    CAMISIM_DIR.mkdir(parents=True, exist_ok=True)

    with open(GENOME_LOCATIONS_PATH, "w") as f:
        for g in genomes:
            fasta_path = (REPO_ROOT / g["local_path"]).resolve()
            f.write(f"{g['genome_ID']}\t{fasta_path}\n")

    with open(METADATA_PATH, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["genome_ID", "OTU", "NCBI_ID", "novelty_category"])
        for g in genomes:
            writer.writerow([g["genome_ID"], g["species_taxid"], g["taxid"], "known_strain"])

    n_otus = len({g["species_taxid"] for g in genomes})
    print(f"Wrote {GENOME_LOCATIONS_PATH} and {METADATA_PATH}")
    print(f"{len(genomes)} genomes across {n_otus} OTUs (species)")


if __name__ == "__main__":
    main()
