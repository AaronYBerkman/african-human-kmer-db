#!/bin/bash
#SBATCH --job-name=per_indiv_afr25
#SBATCH --output=logs/per_indiv_afr25_%a.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=1-25

# Per-individual sensitivity benchmark — African-25 panel.
# 25 continental-African individuals (5 per population: ESN, GWD, LWK, MSL, YRI),
# each with 30,000 reads (15,000 per mate — matched to the Forbes-27 panel's
# depth in 26_per_individual_benchmark.sh, rather than the original 6,000/mate,
# so per-individual estimates are equally precise in both panels). Produces an
# independent Wilcoxon signed-rank test alongside the Forbes-27 result.
#
# Submit from the repository root:
#   sbatch scripts/benchmarks/27_per_individual_african25_benchmark.sh
#
# Outputs: benchmark_results/per_individual_african25_results/<accession>.csv
# Merged:  benchmark_results/per_individual_african25_results/per_individual_african25_summary.csv

source ./config.sh || exit 1   # run from the repository root; sets SCRUBBER, DEFAULT_DB, AFRICAN_DB
READS_DIR="${REPO_ROOT}/reads/human_african25"
META="${REPO_ROOT}/results/human_accessions_african25.tsv"
RESULTS="${REPO_ROOT}/benchmark_results/per_individual_african25_results"
mkdir -p "${RESULTS}"

mapfile -t ACCESSIONS < <(awk -F'\t' 'NR>1 {print $3}' "${META}")
mapfile -t POPULATIONS < <(awk -F'\t' 'NR>1 {print $1}' "${META}")

idx=$((SLURM_ARRAY_TASK_ID - 1))
acc="${ACCESSIONS[$idx]}"
pop="${POPULATIONS[$idx]}"

echo "Task ${SLURM_ARRAY_TASK_ID}: ${pop} ${acc}"

r1="${READS_DIR}/${acc}_1.fastq"
r2="${READS_DIR}/${acc}_2.fastq"

if [ ! -s "${r1}" ] || [ ! -s "${r2}" ]; then
    echo "ERROR: reads missing for ${acc}: ${r1} ${r2}"
    exit 1
fi

# Sequential numeric IDs — ENA mates share the same base read ID, so
# concatenating both files yields duplicate FASTA headers without this.
tmp_fasta=$(mktemp /tmp/per_indiv_afr25_${acc}.XXXXXX.fasta)
awk 'NR%4==1{i++; print ">" i} NR%4==2{print}' "${r1}" "${r2}" > "${tmp_fasta}"

n_reads=$(grep -c "^>" "${tmp_fasta}")
echo "  ${n_reads} reads for ${acc}"

echo "  Running Default..."
${SCRUBBER} aligns_to -db "${DEFAULT_DB}" -num_threads 8 "${tmp_fasta}" \
    2>/dev/null | awk '{print $1}' | sort -u \
    > /tmp/per_indiv_afr25_${acc}_default.txt
default_hits=$(wc -l < /tmp/per_indiv_afr25_${acc}_default.txt)
default_sens=$(awk "BEGIN {printf \"%.6f\", ${default_hits}/${n_reads}}")

echo "  Running African DB..."
${SCRUBBER} aligns_to -db "${AFRICAN_DB}" -num_threads 8 "${tmp_fasta}" \
    2>/dev/null | awk '{print $1}' | sort -u \
    > /tmp/per_indiv_afr25_${acc}_african.txt
african_hits=$(wc -l < /tmp/per_indiv_afr25_${acc}_african.txt)
african_sens=$(awk "BEGIN {printf \"%.6f\", ${african_hits}/${n_reads}}")

# Two-stage combined: Default survivors → African DB
awk 'NR==FNR{hits[$1]=1; next}
     /^>/{id=substr($1,2); skip=(id in hits); if(!skip)print; next}
     !skip{print}' \
    /tmp/per_indiv_afr25_${acc}_default.txt "${tmp_fasta}" \
    > /tmp/per_indiv_afr25_${acc}_survivors.fasta

${SCRUBBER} aligns_to -db "${AFRICAN_DB}" -num_threads 8 /tmp/per_indiv_afr25_${acc}_survivors.fasta \
    2>/dev/null | awk '{print $1}' | sort -u \
    > /tmp/per_indiv_afr25_${acc}_stage2.txt
stage2_hits=$(wc -l < /tmp/per_indiv_afr25_${acc}_stage2.txt)
combined_hits=$((default_hits + stage2_hits))
combined_sens=$(awk "BEGIN {printf \"%.6f\", ${combined_hits}/${n_reads}}")

echo "  Default=${default_sens}  African=${african_sens}  Combined=${combined_sens}"

echo "${pop},${acc},${n_reads},${default_hits},${default_sens},${african_hits},${african_sens},${combined_hits},${combined_sens}" \
    > "${RESULTS}/${acc}.csv"

rm -f "${tmp_fasta}" \
      /tmp/per_indiv_afr25_${acc}_default.txt \
      /tmp/per_indiv_afr25_${acc}_african.txt \
      /tmp/per_indiv_afr25_${acc}_survivors.fasta \
      /tmp/per_indiv_afr25_${acc}_stage2.txt

echo "  Done: ${acc}"

# Task 1 waits for all 25 and merges
if [ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]; then
    echo "Task 1 waiting for all 25 results..."
    while [ "$(ls "${RESULTS}"/*.csv 2>/dev/null | grep -v summary | wc -l)" -lt 25 ]; do
        sleep 20
    done
    echo "population_code,run_accession,n_reads,default_hits,default_sensitivity,african_hits,african_sensitivity,combined_hits,combined_sensitivity" \
        > "${RESULTS}/per_individual_african25_summary.csv"
    cat "${RESULTS}"/ERR*.csv | sort >> "${RESULTS}/per_individual_african25_summary.csv"
    echo "Merged to ${RESULTS}/per_individual_african25_summary.csv"
    cat "${RESULTS}/per_individual_african25_summary.csv"
fi
