#!/bin/bash
#SBATCH --job-name=african25_benchmark
#SBATCH --output=logs/african25_benchmark_%a.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --array=1-18

source ./config.sh || exit 1   # run from the repository root; sets SCRUBBER, DEFAULT_DB, AFRICAN_DB
BENCH_REPO="${REPO_ROOT}/results"
RESULTS="${REPO_ROOT}/benchmark_results/african25_benchmark_results"
mkdir -p $RESULTS

# 18 tasks = 2 arms (microbial, viral) x 9 titres (10-90%)
TASKS=(
    "microbial,${BENCH_REPO}/titrations_microbial_african25/gut_complex_titre10pct_human.fasta,10"
    "microbial,${BENCH_REPO}/titrations_microbial_african25/gut_complex_titre20pct_human.fasta,20"
    "microbial,${BENCH_REPO}/titrations_microbial_african25/gut_complex_titre30pct_human.fasta,30"
    "microbial,${BENCH_REPO}/titrations_microbial_african25/gut_complex_titre40pct_human.fasta,40"
    "microbial,${BENCH_REPO}/titrations_microbial_african25/gut_complex_titre50pct_human.fasta,50"
    "microbial,${BENCH_REPO}/titrations_microbial_african25/gut_complex_titre60pct_human.fasta,60"
    "microbial,${BENCH_REPO}/titrations_microbial_african25/gut_complex_titre70pct_human.fasta,70"
    "microbial,${BENCH_REPO}/titrations_microbial_african25/gut_complex_titre80pct_human.fasta,80"
    "microbial,${BENCH_REPO}/titrations_microbial_african25/gut_complex_titre90pct_human.fasta,90"
    "viral,${BENCH_REPO}/titrations_viral_african25/viral_complex_titre10pct_human.fasta,10"
    "viral,${BENCH_REPO}/titrations_viral_african25/viral_complex_titre20pct_human.fasta,20"
    "viral,${BENCH_REPO}/titrations_viral_african25/viral_complex_titre30pct_human.fasta,30"
    "viral,${BENCH_REPO}/titrations_viral_african25/viral_complex_titre40pct_human.fasta,40"
    "viral,${BENCH_REPO}/titrations_viral_african25/viral_complex_titre50pct_human.fasta,50"
    "viral,${BENCH_REPO}/titrations_viral_african25/viral_complex_titre60pct_human.fasta,60"
    "viral,${BENCH_REPO}/titrations_viral_african25/viral_complex_titre70pct_human.fasta,70"
    "viral,${BENCH_REPO}/titrations_viral_african25/viral_complex_titre80pct_human.fasta,80"
    "viral,${BENCH_REPO}/titrations_viral_african25/viral_complex_titre90pct_human.fasta,90"
)

entry=${TASKS[$SLURM_ARRAY_TASK_ID-1]}
arm=$(echo $entry | cut -d',' -f1)
fasta=$(echo $entry | cut -d',' -f2)
titre=$(echo $entry | cut -d',' -f3)
tag="${arm}_${titre}"

echo "Processing arm=${arm} titre=${titre}% (${fasta})..."

if [ ! -s "$fasta" ]; then
    echo "ERROR: input file missing or empty: $fasta"
    exit 1
fi

tmp_fasta=$(mktemp /tmp/african25_${tag}.XXXXXX.fasta)
cp "$fasta" "$tmp_fasta"

total=$(grep -c "^>" $tmp_fasta)
true_human=$(grep "^>" $tmp_fasta | awk '{print $2}' | grep -c "^ERR")
true_other=$((total - true_human))
human_frac=$(awk "BEGIN {printf \"%.4f\", $true_human/$total}")
echo "  Total: $total | True human: $true_human | True other: $true_other | Human fraction: $human_frac"

# Sort explicitly with -k1,1 and verify before joining anything downstream
# (lesson from the earlier unsorted-join bug: never trust an unverified sort).
grep "^>" $tmp_fasta | awk '$2 ~ /^ERR/{print substr($1,2)}' | sort -k1,1 \
    > /tmp/african25_${tag}_human_ids.txt
echo "  Human ID list built: $(wc -l < /tmp/african25_${tag}_human_ids.txt) IDs"

echo "  Running Condition 1: Default SRA db..."
$SCRUBBER aligns_to -db $DEFAULT_DB -num_threads 8 $tmp_fasta \
    2>/dev/null | awk '{print $1}' | sort -u > /tmp/african25_${tag}_default_hits.txt
default_hits=$(wc -l < /tmp/african25_${tag}_default_hits.txt)
default_tp=$(comm -12 /tmp/african25_${tag}_human_ids.txt /tmp/african25_${tag}_default_hits.txt | wc -l)
default_fp=$((default_hits - default_tp))
default_sens=$(awk "BEGIN {printf \"%.4f\", $default_tp/$true_human}")
default_spec=$(awk "BEGIN {printf \"%.4f\", 1 - $default_fp/$true_other}")
echo "  Default: hits=$default_hits tp=$default_tp fp=$default_fp sens=$default_sens spec=$default_spec"

echo "  Running Condition 2: African db only..."
$SCRUBBER aligns_to -db $AFRICAN_DB -num_threads 8 $tmp_fasta \
    2>/dev/null | awk '{print $1}' | sort -u > /tmp/african25_${tag}_african_hits.txt
african_hits=$(wc -l < /tmp/african25_${tag}_african_hits.txt)
african_tp=$(comm -12 /tmp/african25_${tag}_human_ids.txt /tmp/african25_${tag}_african_hits.txt | wc -l)
african_fp=$((african_hits - african_tp))
african_sens=$(awk "BEGIN {printf \"%.4f\", $african_tp/$true_human}")
african_spec=$(awk "BEGIN {printf \"%.4f\", 1 - $african_fp/$true_other}")
echo "  African: hits=$african_hits tp=$african_tp fp=$african_fp sens=$african_sens spec=$african_spec"

echo "  Running Condition 3: Two-stage pipeline..."
awk 'NR==FNR{hits[$1]=1; next}
     /^>/{
         id=substr($1,2)
         if(id in hits){skip=1} else {skip=0; print}
         next
     }
     !skip{print}' \
    /tmp/african25_${tag}_default_hits.txt $tmp_fasta > /tmp/african25_${tag}_survivors.fasta
echo "  Survivors: $(grep -c "^>" /tmp/african25_${tag}_survivors.fasta) reads"

$SCRUBBER aligns_to -db $AFRICAN_DB -num_threads 8 /tmp/african25_${tag}_survivors.fasta \
    2>/dev/null | awk '{print $1}' | sort -u > /tmp/african25_${tag}_stage2_hits.txt
stage2_hits=$(wc -l < /tmp/african25_${tag}_stage2_hits.txt)
stage2_tp=$(comm -12 /tmp/african25_${tag}_human_ids.txt /tmp/african25_${tag}_stage2_hits.txt | wc -l)
stage2_fp=$((stage2_hits - stage2_tp))
combined_tp=$((default_tp + stage2_tp))
combined_fp=$((default_fp + stage2_fp))
combined_hits=$((default_hits + stage2_hits))
combined_sens=$(awk "BEGIN {printf \"%.4f\", $combined_tp/$true_human}")
combined_spec=$(awk "BEGIN {printf \"%.4f\", 1 - $combined_fp/$true_other}")
echo "  Combined: additional=$stage2_hits total=$combined_hits tp=$combined_tp fp=$combined_fp sens=$combined_sens spec=$combined_spec"

echo "${arm},${titre},$human_frac,$total,$true_human,$true_other,\
$default_hits,$default_tp,$default_fp,$default_sens,$default_spec,\
$african_hits,$african_tp,$african_fp,$african_sens,$african_spec,\
$stage2_hits,$combined_hits,$combined_tp,$combined_fp,$combined_sens,$combined_spec" \
    > $RESULTS/${tag}.csv

rm -f $tmp_fasta /tmp/african25_${tag}_human_ids.txt /tmp/african25_${tag}_default_hits.txt \
      /tmp/african25_${tag}_african_hits.txt /tmp/african25_${tag}_survivors.fasta \
      /tmp/african25_${tag}_stage2_hits.txt

echo "${tag} complete."

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
    echo "Task 1 waiting for all 18 results..."
    while [ $(ls $RESULTS/*.csv 2>/dev/null | grep -E "(microbial|viral)_[0-9]+\.csv" | wc -l) -lt 18 ]; do
        sleep 30
    done
    echo "arm,human_titre_pct,human_fraction,total_reads,true_human,true_other,\
default_hits,default_tp,default_fp,default_sensitivity,default_specificity,\
african_hits,african_tp,african_fp,african_sensitivity,african_specificity,\
stage2_additional_hits,combined_hits,combined_tp,combined_fp,combined_sensitivity,combined_specificity" \
        > $RESULTS/african25_benchmark_summary.csv
    sort -t, -k1,1 -k2,2n $RESULTS/microbial_*.csv $RESULTS/viral_*.csv >> $RESULTS/african25_benchmark_summary.csv
    echo "Merge complete."
    cat $RESULTS/african25_benchmark_summary.csv
fi
