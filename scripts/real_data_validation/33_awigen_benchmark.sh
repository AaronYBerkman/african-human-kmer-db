#!/bin/bash
#SBATCH --job-name=awigen_benchmark
#SBATCH --output=logs/awigen_benchmark_%a.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-600%50

# The container, database, raw reads and outputs must be on storage shared
# by all compute nodes. Check a few tasks' outputs before trusting a full
# run: in our runs, one node returned zero hits for both databases on every
# task scheduled there, and those tasks had to be rerun elsewhere.

# Real-data validation: does the African database catch more human reads than
# the default on genuine AWI-Gen 2 metagenomic sequencing? Unlike the
# synthetic titration benchmark, there is no per-read ground truth here (real
# sequencer headers don't say "this read is human") - this measures
# incremental capture (how many more reads each condition flags as human),
# not sensitivity/specificity.
#
# Two conditions only - Default and African DB, each run standalone on the
# full read set - not the three-condition (+Combined two-stage) design used
# in the synthetic benchmark. The comparison here is simply how many reads
# each database removes on its own.
#
# Full depth, no subsampling: every read from every lane for each selected
# individual. Each raw file (every lane, both mates) is classified
# separately and hit counts are summed - NOT concatenated into one FASTA
# first. Concatenating multiple gzip-sourced files into a single FASTA
# before calling aligns_to was tried first and produces silently wrong
# results (zero hits from both Default and African DB on a real 44.5M-read,
# 4-file individual, despite aligns_to finding real hits on single-file
# tests up to 5M reads) - root cause not fully diagnosed, but isolated
# specifically to multi-file concatenation, not scale (a 5M-read single-file
# test worked fine; a 4-file test failed even at just 4K reads). Per-file
# processing sidesteps the bug entirely: each raw file is a physically
# distinct set of reads (different lane/mate), so summing independently
# deduplicated per-file hit counts gives the correct total with no need to
# merge read IDs across files.
#
# Submit from the repository root (needs results/awigen_validation_samples.tsv
# from 32_select_awigen_samples.py, and AWIGEN_RAW_DIR set in config.sh):
#   sbatch scripts/real_data_validation/33_awigen_benchmark.sh
#
# Outputs: benchmark_results/awigen_benchmark_results/<metaG_id>.csv, one per
# individual (not distributed). No merge step: 34_awigen_analysis.R reads the
# per-individual files directly once all 600 tasks are confirmed complete.

source ./config.sh || exit 1   # run from the repository root; sets SCRUBBER, DEFAULT_DB, AFRICAN_DB, AWIGEN_RAW_DIR
RAW_DATA_DIR="${AWIGEN_RAW_DIR}"
MANIFEST="${REPO_ROOT}/results/awigen_validation_samples.tsv"
RESULTS="${REPO_ROOT}/benchmark_results/awigen_benchmark_results"
mkdir -p "${RESULTS}"

# Manifest has a header row, so data row N is line N+1
metaG_id=$(awk -F'\t' -v idx=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==idx {print $1}' "${MANIFEST}")
if [ -z "${metaG_id}" ]; then
    echo "ERROR: no metaG_id found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi
echo "Task ${SLURM_ARRAY_TASK_ID}: ${metaG_id}"

if [ -s "${RESULTS}/${metaG_id}.csv" ]; then
    echo "${metaG_id}.csv already exists, skipping (resume-safe)."
    exit 0
fi

all_files=(${RAW_DATA_DIR}/${metaG_id}_*_1.fq.gz ${RAW_DATA_DIR}/${metaG_id}_*_2.fq.gz)
if [ ! -s "${all_files[0]}" ]; then
    echo "ERROR: no raw data files found for ${metaG_id} in ${RAW_DATA_DIR}"
    exit 1
fi
echo "  ${#all_files[@]} file(s): ${all_files[*]##*/}"

total=0
default_hits=0
african_hits=0

for f in "${all_files[@]}"; do
    fname=$(basename "$f")
    tmp_fasta=$(mktemp /tmp/awigen_${metaG_id}.XXXXXX.fasta)
    zcat "$f" | awk 'NR%4==1{i++; print ">" i} NR%4==2{print}' > "${tmp_fasta}"

    n_reads=$(grep -c "^>" "${tmp_fasta}")

    n_default=$(${SCRUBBER} aligns_to -db "${DEFAULT_DB}" -num_threads 8 "${tmp_fasta}" \
        2>/dev/null | awk '{print $1}' | sort -u | wc -l)

    n_african=$(${SCRUBBER} aligns_to -db "${AFRICAN_DB}" -num_threads 8 "${tmp_fasta}" \
        2>/dev/null | awk '{print $1}' | sort -u | wc -l)

    echo "  ${fname}: ${n_reads} reads, default=${n_default}, african=${n_african}"

    total=$((total + n_reads))
    default_hits=$((default_hits + n_default))
    african_hits=$((african_hits + n_african))

    rm -f "${tmp_fasta}"
done

echo "  TOTAL: ${total} reads  Default=${default_hits}  African=${african_hits}"
echo "  African removes $((african_hits - default_hits)) more reads than Default."

echo "${metaG_id},${total},${default_hits},${african_hits}" \
    > "${RESULTS}/${metaG_id}.csv"

echo "${metaG_id} complete."
