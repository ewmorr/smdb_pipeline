#!/bin/bash
# =====================================================================
# Submit the whole pipeline with dependencies:
#   step 1 (Kraken2, big-memory array)
#     -> step 2 (Architeuthis + Bracken, small array; task N starts when Kraken2 task N succeeds)
#       -> step 3 (summary, once all step-2 tasks succeed)
#
#   bash submit_all.sh              # all samples
#   bash submit_all.sh 1-3          # only samples 1-3 of samples.tsv (good first test)
#   bash submit_all.sh 4,7          # specific samples
# =====================================================================
set -euo pipefail
PIPE_DIR="$(cd "$(dirname "$0")" && pwd)"
export PIPE_DIR
source "$PIPE_DIR/config.sh"

[[ -s "$SAMPLES_TSV" ]] || { echo "Run 'bash 00_setup.sh' first."; exit 1; }
[[ -x "$KRAKEN2_REPORT_BIN" ]] || { echo "kraken2-report missing; run 'bash 00_setup.sh' first."; exit 1; }
NSAMP=$(( $(wc -l < "$SAMPLES_TSV") - 1 ))
RANGE="${1:-1-$NSAMP}"
mkdir -p "$LOG_DIR"

# Optional sbatch flags built from config.sh
k_opts=(); f_opts=()
[[ -n "$KRAKEN_PARTITION"  ]] && k_opts+=(--partition="$KRAKEN_PARTITION")
[[ -n "$KRAKEN_CONSTRAINT" ]] && k_opts+=(--constraint="$KRAKEN_CONSTRAINT")
[[ -n "$FILTER_PARTITION"  ]] && f_opts+=(--partition="$FILTER_PARTITION")
# shellcheck disable=SC2206
common=($SBATCH_COMMON)

echo "Samples: $RANGE of $NSAMP   |   Kraken2 mode: $KRAKEN_MODE, ${KRAKEN_MEM}, ${KRAKEN_CPUS} CPUs"

J1=$(sbatch --parsable ${common[@]+"${common[@]}"} ${k_opts[@]+"${k_opts[@]}"} \
    --array="${RANGE}%${KRAKEN_MAX_PARALLEL}" \
    --mem="$KRAKEN_MEM" --cpus-per-task="$KRAKEN_CPUS" --time="$KRAKEN_TIME" \
    --output="$LOG_DIR/kraken2_%A_%a.log" \
    --export=ALL,PIPE_DIR="$PIPE_DIR" \
    "$PIPE_DIR/01_kraken2.slurm")
echo "Step 1 Kraken2:              job $J1"

J2=$(sbatch --parsable ${common[@]+"${common[@]}"} ${f_opts[@]+"${f_opts[@]}"} \
    --array="$RANGE" \
    --dependency="aftercorr:$J1" \
    --mem="$FILTER_MEM" --cpus-per-task="$FILTER_CPUS" --time="$FILTER_TIME" \
    --output="$LOG_DIR/filter_bracken_%A_%a.log" \
    --export=ALL,PIPE_DIR="$PIPE_DIR" \
    "$PIPE_DIR/02_filter_bracken.slurm")
echo "Step 2 Architeuthis+Bracken: job $J2 (each task waits for its Kraken2 task)"

J3=$(sbatch --parsable ${common[@]+"${common[@]}"} ${f_opts[@]+"${f_opts[@]}"} \
    --dependency="afterok:$J2" \
    --output="$LOG_DIR/summary_%j.log" \
    --export=ALL,PIPE_DIR="$PIPE_DIR" \
    "$PIPE_DIR/03_summarize.slurm")
echo "Step 3 summary:              job $J3"
echo
echo "Watch with:  squeue -u $USER      Logs in: $LOG_DIR"
