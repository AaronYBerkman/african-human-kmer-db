#!/bin/bash
#SBATCH --job-name=decompress
#SBATCH --output=logs/decompress.log
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00

source ./config.sh || exit 1   # run from the repository root

cd ${WORK_DIR}/microbial_genomes/

for f in *.fna.gz; do
    gunzip -k "$f"
done

echo "Decompression complete"
