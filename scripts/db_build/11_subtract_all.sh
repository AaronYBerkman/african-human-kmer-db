#!/bin/bash
#SBATCH --job-name=subtract_all
#SBATCH --output=logs/subtract_all.log
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=12:00:00

source ./config.sh || exit 1   # run from the repository root

# Step 1: Subtract microbial k-mers from African human database
${STAT_BIN}/subtract_db \
  ${STAT_DB_DIR}/african_human.dense.db \
  ${STAT_DB_DIR}/microbial.dense.db \
  ${STAT_DB_DIR}/african_human_minus_microbial.db

# Step 2: Subtract viral k-mers from result of Step 1
${STAT_BIN}/subtract_db \
  ${STAT_DB_DIR}/african_human_minus_microbial.db \
  ${STAT_DB_DIR}/viral.dense.db \
  ${STAT_DB_DIR}/african_human_filtered.db
