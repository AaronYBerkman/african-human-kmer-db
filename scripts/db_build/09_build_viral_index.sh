#!/bin/bash
#SBATCH --job-name=build_viral_index
#SBATCH --output=logs/build_viral_index_%a.log
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=168:00:00
#SBATCH --array=1-16

source ./config.sh || exit 1   # run from the repository root
# viral_split_lists/viral_files_*: the viral genome list split into 16 chunks,
# one per array task (written by hand; see README, Database construction).

CHUNK=$(ls ${STAT_DB_DIR}/viral_split_lists/viral_files_* | sed -n "${SLURM_ARRAY_TASK_ID}p")

${STAT_BIN}/build_index_of_each_file \
  $CHUNK \
  4 \
  32 \
  .dense.db
