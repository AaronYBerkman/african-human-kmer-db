#!/bin/bash
#SBATCH --job-name=merge_microbial
#SBATCH --output=logs/merge_microbial.log
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00

source ./config.sh || exit 1   # run from the repository root
# microbial_files.list: file size <tab> path for every microbial .fna, the same
# format 10_merge_viral_db.sh writes for the viral genomes.

${STAT_BIN}/merge_db \
  ${STAT_DB_DIR}/microbial_files.list \
  .dense.db \
  ${STAT_DB_DIR}/microbial.dense.db \
  30000000000
