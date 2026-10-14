#!/bin/bash
#SBATCH --job-name=viral_download
#SBATCH --output=logs/viral_download.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=168:00:00

source ./config.sh || exit 1   # run from the repository root

mkdir -p ${WORK_DIR}/viral_genomes

awk -F'\t' '$12=="Complete Genome" && $11=="latest" {print $20}' \
  ${WORK_DIR}/taxonomy/viral_assembly_summary.txt | \
  while read url; do
    fname=$(basename $url)
    wget -q -c "${url}/${fname}_genomic.fna.gz" \
      -O "${WORK_DIR}/viral_genomes/${fname}_genomic.fna.gz"
  done

echo "Download complete"
