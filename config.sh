# Paths for the steps that run on a compute cluster: database construction,
# benchmark runs and the AWI-Gen 2 validation. Edit these for your system and
# run every script from the repository root, e.g.
#   sbatch scripts/db_build/01_build_index.sh
# The analysis scripts (scripts/analysis/) do not use this file.

REPO_ROOT="${REPO_ROOT:-$(pwd)}"

# Large working data: downloaded genomes, k-mer indexes and the built database
WORK_DIR="${WORK_DIR:-${REPO_ROOT}/work}"
STAT_DB_DIR="${STAT_DB_DIR:-${WORK_DIR}/stat_db}"

# NCBI STAT (ngs-tools, branch "tax"): directory holding build_index_of_each_file,
# merge_db and subtract_db
STAT_BIN="${STAT_BIN:-/path/to/ngs-tools/tools/tax/bin}"

# NCBI SRA Human Scrubber 2.2.1 container and the two databases compared
SCRUBBER_SIF="${SCRUBBER_SIF:-/path/to/human_scrubber_2.2.1.sif}"
SCRUBBER="singularity exec --bind ${REPO_ROOT},${WORK_DIR} ${SCRUBBER_SIF}"
DEFAULT_DB="/opt/scrubber/data/human_filter.db"   # default database, inside the container
AFRICAN_DB="${AFRICAN_DB:-${STAT_DB_DIR}/african_human_filtered.db}"

# AWI-Gen 2 validation only: raw sequencing reads (controlled access, not distributed)
AWIGEN_RAW_DIR="${AWIGEN_RAW_DIR:-/path/to/awigen2_raw_reads}"
