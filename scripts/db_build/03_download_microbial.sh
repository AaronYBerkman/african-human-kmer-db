#!/bin/bash
#SBATCH --job-name=microbial_download
#SBATCH --output=logs/microbial_download.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00

source ./config.sh || exit 1   # run from the repository root

mkdir -p ${WORK_DIR}/microbial_genomes

# Download bacterial complete genomes
awk -F'\t' '$12=="Complete Genome" && $11=="latest" {print $20}' \
  ${WORK_DIR}/taxonomy/bacterial_assembly_summary.txt | \
  while read url; do
    fname=$(basename $url)
    wget -q -c "${url}/${fname}_genomic.fna.gz" \
      -O "${WORK_DIR}/microbial_genomes/${fname}_genomic.fna.gz"
  done

# Download archaeal complete genomes
awk -F'\t' '$12=="Complete Genome" && $11=="latest" {print $20}' \
  ${WORK_DIR}/taxonomy/archaeal_assembly_summary.txt | \
  while read url; do
    fname=$(basename $url)
    wget -q -c "${url}/${fname}_genomic.fna.gz" \
      -O "${WORK_DIR}/microbial_genomes/${fname}_genomic.fna.gz"
  done

echo "Download complete"
