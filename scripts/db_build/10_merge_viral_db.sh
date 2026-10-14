#!/bin/bash
#SBATCH --job-name=merge_viral_db
#SBATCH --output=logs/merge_viral_db.log
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=6:00:00

source ./config.sh || exit 1   # run from the repository root

# Generate files.list: filesize + path to each .fna (merge_db appends .dense.db)
find ${WORK_DIR}/viral_genomes/ -name "*.fna" | while read f; do
    size=$(stat -c '%s' "${f}.dense.db")
    echo -e "${size}\t${f}"
done > ${STAT_DB_DIR}/viral_files.list

# Merge all viral .dense.db files into one
${STAT_BIN}/merge_db \
  ${STAT_DB_DIR}/viral_files.list \
  .dense.db \
  ${STAT_DB_DIR}/viral.dense.db \
  300000000
