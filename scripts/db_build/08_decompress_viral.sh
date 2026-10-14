#!/bin/bash
#SBATCH --job-name=decompress_viral
#SBATCH --output=logs/decompress_viral.log
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=4:00:00

source ./config.sh || exit 1   # run from the repository root

ls ${WORK_DIR}/viral_genomes/*.fna.gz | xargs -P 16 -I {} gunzip -k {}

echo "Decompression complete"
