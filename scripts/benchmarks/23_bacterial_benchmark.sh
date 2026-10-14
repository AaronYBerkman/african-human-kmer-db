#!/bin/bash
#SBATCH --job-name=bacterial_benchmark
#SBATCH --output=logs/bacterial_benchmark_%a.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --array=1-9

source ./config.sh || exit 1   # run from the repository root; sets SCRUBBER, DEFAULT_DB, AFRICAN_DB
INPUT_DIR="${REPO_ROOT}/results/titrations"
RESULTS="${REPO_ROOT}/benchmark_results/bacterial_benchmark_results"
mkdir -p $RESULTS

# 9 titration points, 10-90% human, complex 29-genome bacterial background
TITRES=(10 20 30 40 50 60 70 80 90)

titre=${TITRES[$SLURM_ARRAY_TASK_ID-1]}
fasta="${INPUT_DIR}/gut_complex_titre${titre}pct_human.fasta"

echo "Processing titre=${titre}% (${fasta})..."

if [ ! -s "$fasta" ]; then
    echo "ERROR: input file missing or empty: $fasta"
    exit 1
fi

# Already FASTA (not FASTQ) - no conversion needed, but copy to scratch tmp
# so we don't repeatedly read the shared filesystem for every grep/awk pass below.
tmp_fasta=$(mktemp /tmp/gut_complex_${titre}pct.XXXXXX.fasta)
cp "$fasta" "$tmp_fasta"

# Ground truth from read names: field 2 starting with "ERR" = human,
# anything else (RefSeq genome accession, e.g. GCF_...) = bacterial.
# Same tagging convention used by every benchmark script in this repo.
total=$(grep -c "^>" $tmp_fasta)
true_human=$(grep "^>" $tmp_fasta | awk '{print $2}' | grep -c "^ERR")
true_bacterial=$((total - true_human))
human_frac=$(awk "BEGIN {printf \"%.4f\", $true_human/$total}")
echo "  Total: $total | True human: $true_human | True bacterial: $true_bacterial | Human fraction: $human_frac"

# Build sorted human read ID list for TP calculation
grep "^>" $tmp_fasta | awk '$2 ~ /^ERR/{print substr($1,2)}' | sort \
    > /tmp/gut_complex_${titre}pct_human_ids.txt
echo "  Human ID list built: $(wc -l < /tmp/gut_complex_${titre}pct_human_ids.txt) IDs"

# Condition 1: Default SRA db only
echo "  Running Condition 1: Default SRA db..."
$SCRUBBER aligns_to -db $DEFAULT_DB -num_threads 8 $tmp_fasta \
    2>/dev/null | awk '{print $1}' | sort -u > /tmp/gut_complex_${titre}pct_default_hits.txt
default_hits=$(wc -l < /tmp/gut_complex_${titre}pct_default_hits.txt)
default_tp=$(comm -12 /tmp/gut_complex_${titre}pct_human_ids.txt /tmp/gut_complex_${titre}pct_default_hits.txt | wc -l)
default_fp=$((default_hits - default_tp))
default_sens=$(awk "BEGIN {printf \"%.4f\", $default_tp/$true_human}")
default_spec=$(awk "BEGIN {printf \"%.4f\", 1 - $default_fp/$true_bacterial}")
echo "  Default: hits=$default_hits tp=$default_tp fp=$default_fp sens=$default_sens spec=$default_spec"

# Condition 2: African db only
echo "  Running Condition 2: African db only..."
$SCRUBBER aligns_to -db $AFRICAN_DB -num_threads 8 $tmp_fasta \
    2>/dev/null | awk '{print $1}' | sort -u > /tmp/gut_complex_${titre}pct_african_hits.txt
african_hits=$(wc -l < /tmp/gut_complex_${titre}pct_african_hits.txt)
african_tp=$(comm -12 /tmp/gut_complex_${titre}pct_human_ids.txt /tmp/gut_complex_${titre}pct_african_hits.txt | wc -l)
african_fp=$((african_hits - african_tp))
african_sens=$(awk "BEGIN {printf \"%.4f\", $african_tp/$true_human}")
african_spec=$(awk "BEGIN {printf \"%.4f\", 1 - $african_fp/$true_bacterial}")
echo "  African: hits=$african_hits tp=$african_tp fp=$african_fp sens=$african_sens spec=$african_spec"

# Condition 3: Two-stage
echo "  Running Condition 3: Two-stage pipeline..."
awk 'NR==FNR{hits[$1]=1; next}
     /^>/{
         id=substr($1,2)
         if(id in hits){skip=1} else {skip=0; print}
         next
     }
     !skip{print}' \
    /tmp/gut_complex_${titre}pct_default_hits.txt $tmp_fasta > /tmp/gut_complex_${titre}pct_survivors.fasta
echo "  Survivors: $(grep -c "^>" /tmp/gut_complex_${titre}pct_survivors.fasta) reads"

$SCRUBBER aligns_to -db $AFRICAN_DB -num_threads 8 /tmp/gut_complex_${titre}pct_survivors.fasta \
    2>/dev/null | awk '{print $1}' | sort -u > /tmp/gut_complex_${titre}pct_stage2_hits.txt
stage2_hits=$(wc -l < /tmp/gut_complex_${titre}pct_stage2_hits.txt)
stage2_tp=$(comm -12 /tmp/gut_complex_${titre}pct_human_ids.txt /tmp/gut_complex_${titre}pct_stage2_hits.txt | wc -l)
stage2_fp=$((stage2_hits - stage2_tp))
combined_tp=$((default_tp + stage2_tp))
combined_fp=$((default_fp + stage2_fp))
combined_hits=$((default_hits + stage2_hits))
combined_sens=$(awk "BEGIN {printf \"%.4f\", $combined_tp/$true_human}")
combined_spec=$(awk "BEGIN {printf \"%.4f\", 1 - $combined_fp/$true_bacterial}")
echo "  Combined: additional=$stage2_hits total=$combined_hits tp=$combined_tp fp=$combined_fp sens=$combined_sens spec=$combined_spec"

# Write CSV row
echo "${titre},$human_frac,$total,$true_human,$true_bacterial,\
$default_hits,$default_tp,$default_fp,$default_sens,$default_spec,\
$african_hits,$african_tp,$african_fp,$african_sens,$african_spec,\
$stage2_hits,$combined_hits,$combined_tp,$combined_fp,$combined_sens,$combined_spec" \
    > $RESULTS/titre${titre}.csv

# Cleanup
rm -f $tmp_fasta /tmp/gut_complex_${titre}pct_human_ids.txt /tmp/gut_complex_${titre}pct_default_hits.txt \
      /tmp/gut_complex_${titre}pct_african_hits.txt /tmp/gut_complex_${titre}pct_survivors.fasta \
      /tmp/gut_complex_${titre}pct_stage2_hits.txt

echo "titre=${titre}% complete."

# Task 1 waits for all 9 results then merges
if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
    echo "Task 1 waiting for all 9 results..."
    while [ $(ls $RESULTS/titre*.csv 2>/dev/null | wc -l) -lt 9 ]; do
        sleep 30
    done
    echo "human_titre_pct,human_fraction,total_reads,true_human,true_bacterial,\
default_hits,default_tp,default_fp,default_sensitivity,default_specificity,\
african_hits,african_tp,african_fp,african_sensitivity,african_specificity,\
stage2_additional_hits,combined_hits,combined_tp,combined_fp,combined_sensitivity,combined_specificity" \
        > $RESULTS/bacterial_benchmark_summary.csv
    sort -t, -k1,1n $RESULTS/titre*.csv >> $RESULTS/bacterial_benchmark_summary.csv
    echo "Merge complete."
    cat $RESULTS/bacterial_benchmark_summary.csv
fi
