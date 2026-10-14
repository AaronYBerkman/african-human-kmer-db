#!/bin/bash
#SBATCH --job-name=merge_db
#SBATCH --output=logs/merge_db.log
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --time=6:00:00

source ./config.sh || exit 1   # run from the repository root

cd ${STAT_DB_DIR} && \
  ${STAT_BIN}/merge_db \
  ./files.list \
  .dense.db \
  ./african_human.dense.db \
  9000000000
