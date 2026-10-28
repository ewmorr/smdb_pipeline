#!/bin/bash
# =====================================================================
# Step 1: Kraken2 classification against SoilMicrobeDB (one array task per sample)
# Submitted by submit_all.sh, which sets memory/CPU/partition from config.sh.
# Manual use:  sbatch --array=1-19%3 --mem=340G --cpus-per-task=16 01_kraken2.slurm
# =====================================================================
#SBATCH --job-name=smdb_kraken2
#SBATCH --output=%x_%A_%a.log
set -euo pipefail

source "${PIPE_DIR:-$(dirname "$0")}/config.sh"
conda_init

load_sample "$SLURM_ARRAY_TASK_ID"
THREADS="${SLURM_CPUS_PER_TASK:-$KRAKEN_CPUS}"
REPORT="$KRAKEN_OUT_DIR/${SAMPLE}_kraken.kreport"
OUTPUT="$KRAKEN_OUT_DIR/${SAMPLE}_kraken.output"

echo "[$(date)] Sample $SAMPLE on $(hostname), $THREADS threads, mode=$KRAKEN_MODE"
echo "Input: $FASTQ"

if [[ -s "$REPORT" && -s "$OUTPUT" && -f "$KRAKEN_OUT_DIR/${SAMPLE}.done" ]]; then
    echo "Already finished; skipping."
    exit 0
fi

WORK="$(mktemp -d "$LOCAL_SCRATCH/smdb_${SAMPLE}_XXXX")"
trap 'rm -rf "$WORK"' EXIT

# --- JGI files are interleaved; Kraken2 --paired needs separate mate files
echo "[$(date)] De-interleaving into $WORK"
reformat.sh -Xmx4g in="$FASTQ" out1="$WORK/R1.fq.gz" out2="$WORK/R2.fq.gz" \
    int=t verifypaired=t zl=2 ow=t

# --- Choose how the ~280 GB hash is accessed
DB_USE="$DBDIR"
MMAP_FLAG=()
if [[ "$KRAKEN_MODE" == "mmap" ]]; then
    MMAP_FLAG=(--memory-mapping)
    if [[ "$MMAP_COPY_TO_LOCAL" == "true" ]]; then
        echo "[$(date)] Copying database to node-local disk"
        mkdir -p "$WORK/db"
        cp "$DBDIR"/{hash.k2d,opts.k2d,taxo.k2d} "$WORK/db/"
        DB_USE="$WORK/db"
    fi
fi

echo "[$(date)] Running Kraken2"
kraken2 --db "$DB_USE" \
    --threads "$THREADS" \
    ${MMAP_FLAG[@]+"${MMAP_FLAG[@]}"} \
    --report "$REPORT.tmp" \
    --output "$OUTPUT.tmp" \
    --report-zero-counts --report-minimizer-data \
    --paired "$WORK/R1.fq.gz" "$WORK/R2.fq.gz"

mv "$REPORT.tmp" "$REPORT"
mv "$OUTPUT.tmp" "$OUTPUT"
touch "$KRAKEN_OUT_DIR/${SAMPLE}.done"
echo "[$(date)] Done: $REPORT"
