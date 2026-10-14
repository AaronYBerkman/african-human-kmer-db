#!/bin/bash
#SBATCH --job-name=microbial_index
#SBATCH --output=logs/microbial_index_%a.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=168:00:00
#SBATCH --array=1-16

source ./config.sh || exit 1   # run from the repository root
# split_lists/microbial_files_*: the microbial genome list split into 16 chunks,
# one per array task (written by hand; see README, Database construction).

# Get the chunk file for this array task
CHUNK=$(ls ${STAT_DB_DIR}/split_lists/microbial_files_* | sed -n "${SLURM_ARRAY_TASK_ID}p")

${STAT_BIN}/build_index_of_each_file \
  $CHUNK \
  4 \
  32 \
  .dense.db
