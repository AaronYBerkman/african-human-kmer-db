#!/bin/bash
#SBATCH --job-name=build_index
#SBATCH --output=logs/build_index.log
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00

source ./config.sh || exit 1   # run from the repository root
# files.list (file size <tab> path, one line per haplotype FASTA of the three
# African assemblies) is written by hand; see README, Database construction.

cd ${STAT_DB_DIR} && \
  ${STAT_BIN}/build_index_of_each_file \
  ./files.list \
  4 \
  32 \
  .dense.db
