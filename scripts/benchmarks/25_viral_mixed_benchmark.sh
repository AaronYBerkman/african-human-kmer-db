#!/bin/bash
#SBATCH --job-name=viral_mixed_benchmark
#SBATCH --output=logs/viral_mixed_benchmark_%a.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --array=1-9

# Viral x Mixed-ancestry benchmark: 200,000 reads per file, head-of-file
# human reads from the 27 Forbes-27 accessions (human_pool.fasta), pure
# viral reads from the 27 ERP169281 source accessions (viral_pool.fasta,
# built by 20_extract_viral_pool.py). Human reads are drawn the same way as
# the microbial arm rather than using Forbes' own pre-mixed files directly,
# which embed a human background of different, uncontrolled composition.
#
# Prerequisite (from the repository root):
#   python scripts/titration_gen/17_generate_titrations.py \
#     --content-pool results/pools/viral_pool.fasta \
#     --human-pool   results/pools/human_pool.fasta \
#     --output-dir   results/titrations_viral_forbes27 \
#     --filename-prefix viral_complex \
#     --manifest-path results/manifest_viral_forbes27.tsv \
#     --content-label viral
#
# Then submit from the repository root:
#   sbatch scripts/benchmarks/25_viral_mixed_benchmark.sh

source ./config.sh || exit 1   # run from the repository root; sets SCRUBBER, DEFAULT_DB, AFRICAN_DB
BENCH_REPO="${REPO_ROOT}/results"
RESULTS="${REPO_ROOT}/benchmark_results/viral_mixed_benchmark_results"
mkdir -p "${RESULTS}"

TITRES=(10 20 30 40 50 60 70 80 90)
titre=${TITRES[$SLURM_ARRAY_TASK_ID-1]}
fasta="${BENCH_REPO}/titrations_viral_forbes27/viral_complex_titre${titre}pct_human.fasta"
tag="viral_${titre}"

echo "Processing titre=${titre}% (${fasta})..."

if [ ! -s "${fasta}" ]; then
    echo "ERROR: input file missing or empty: ${fasta}"
    echo "Run scripts/titration_gen/17_generate_titrations.py first."
    exit 1
fi

# Copy to local /tmp to avoid repeated NFS reads across the grep/awk passes below
tmp_fasta=$(mktemp /tmp/viral_mixed_${tag}.XXXXXX.fasta)
cp "${fasta}" "${tmp_fasta}"

total=$(grep -c "^>" "${tmp_fasta}")
true_human=$(grep "^>" "${tmp_fasta}" | awk '{print $2}' | grep -c "^ERR")
true_other=$((total - true_human))
human_frac=$(awk "BEGIN {printf \"%.4f\", ${true_human}/${total}}")
echo "  Total: ${total} | True human: ${true_human} | True viral: ${true_other} | Human fraction: ${human_frac}"

grep "^>" "${tmp_fasta}" | awk '$2 ~ /^ERR/{print substr($1,2)}' | sort -k1,1 \
    > /tmp/viral_mixed_${tag}_human_ids.txt
echo "  Human ID list: $(wc -l < /tmp/viral_mixed_${tag}_human_ids.txt) IDs"

echo "  Running Condition 1: Default..."
${SCRUBBER} aligns_to -db "${DEFAULT_DB}" -num_threads 8 "${tmp_fasta}" \
    2>/dev/null | awk '{print $1}' | sort -u > /tmp/viral_mixed_${tag}_default_hits.txt
default_hits=$(wc -l < /tmp/viral_mixed_${tag}_default_hits.txt)
default_tp=$(comm -12 /tmp/viral_mixed_${tag}_human_ids.txt /tmp/viral_mixed_${tag}_default_hits.txt | wc -l)
default_fp=$((default_hits - default_tp))
default_sens=$(awk "BEGIN {printf \"%.4f\", ${default_tp}/${true_human}}")
default_spec=$(awk "BEGIN {printf \"%.4f\", 1 - ${default_fp}/${true_other}}")
echo "  Default: hits=${default_hits} tp=${default_tp} fp=${default_fp} sens=${default_sens} spec=${default_spec}"

echo "  Running Condition 2: African DB only..."
${SCRUBBER} aligns_to -db "${AFRICAN_DB}" -num_threads 8 "${tmp_fasta}" \
    2>/dev/null | awk '{print $1}' | sort -u > /tmp/viral_mixed_${tag}_african_hits.txt
african_hits=$(wc -l < /tmp/viral_mixed_${tag}_african_hits.txt)
african_tp=$(comm -12 /tmp/viral_mixed_${tag}_human_ids.txt /tmp/viral_mixed_${tag}_african_hits.txt | wc -l)
african_fp=$((african_hits - african_tp))
african_sens=$(awk "BEGIN {printf \"%.4f\", ${african_tp}/${true_human}}")
african_spec=$(awk "BEGIN {printf \"%.4f\", 1 - ${african_fp}/${true_other}}")
echo "  African: hits=${african_hits} tp=${african_tp} fp=${african_fp} sens=${african_sens} spec=${african_spec}"

echo "  Running Condition 3: Two-stage pipeline..."
awk 'NR==FNR{hits[$1]=1; next}
     /^>/{id=substr($1,2); if(id in hits){skip=1} else {skip=0; print}; next}
     !skip{print}' \
    /tmp/viral_mixed_${tag}_default_hits.txt "${tmp_fasta}" \
    > /tmp/viral_mixed_${tag}_survivors.fasta

${SCRUBBER} aligns_to -db "${AFRICAN_DB}" -num_threads 8 /tmp/viral_mixed_${tag}_survivors.fasta \
    2>/dev/null | awk '{print $1}' | sort -u > /tmp/viral_mixed_${tag}_stage2_hits.txt
stage2_hits=$(wc -l < /tmp/viral_mixed_${tag}_stage2_hits.txt)
stage2_tp=$(comm -12 /tmp/viral_mixed_${tag}_human_ids.txt /tmp/viral_mixed_${tag}_stage2_hits.txt | wc -l)
stage2_fp=$((stage2_hits - stage2_tp))
combined_tp=$((default_tp + stage2_tp))
combined_fp=$((default_fp + stage2_fp))
combined_hits=$((default_hits + stage2_hits))
combined_sens=$(awk "BEGIN {printf \"%.4f\", ${combined_tp}/${true_human}}")
combined_spec=$(awk "BEGIN {printf \"%.4f\", 1 - ${combined_fp}/${true_other}}")
echo "  Combined: additional=${stage2_hits} tp=${combined_tp} fp=${combined_fp} sens=${combined_sens} spec=${combined_spec}"

echo "viral,${titre},${human_frac},${total},${true_human},${true_other},\
${default_hits},${default_tp},${default_fp},${default_sens},${default_spec},\
${african_hits},${african_tp},${african_fp},${african_sens},${african_spec},\
${stage2_hits},${combined_hits},${combined_tp},${combined_fp},${combined_sens},${combined_spec}" \
    > "${RESULTS}/${tag}.csv"

rm -f "${tmp_fasta}" \
      /tmp/viral_mixed_${tag}_human_ids.txt \
      /tmp/viral_mixed_${tag}_default_hits.txt \
      /tmp/viral_mixed_${tag}_african_hits.txt \
      /tmp/viral_mixed_${tag}_survivors.fasta \
      /tmp/viral_mixed_${tag}_stage2_hits.txt

echo "${tag} complete."

# Task 1 waits for all 9 and merges
if [ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]; then
    echo "Task 1 waiting for all 9 results..."
    while [ "$(ls "${RESULTS}"/viral_*.csv 2>/dev/null | wc -l)" -lt 9 ]; do
        sleep 30
    done
    echo "arm,human_titre_pct,human_fraction,total_reads,true_human,true_other,\
default_hits,default_tp,default_fp,default_sensitivity,default_specificity,\
african_hits,african_tp,african_fp,african_sensitivity,african_specificity,\
stage2_additional_hits,combined_hits,combined_tp,combined_fp,combined_sensitivity,combined_specificity" \
        > "${RESULTS}/viral_mixed_benchmark_summary.csv"
    sort -t, -k2,2n "${RESULTS}"/viral_[0-9]*.csv >> "${RESULTS}/viral_mixed_benchmark_summary.csv"
    echo "Merged to ${RESULTS}/viral_mixed_benchmark_summary.csv"
    cat "${RESULTS}/viral_mixed_benchmark_summary.csv"
fi
