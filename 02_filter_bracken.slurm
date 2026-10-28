#!/bin/bash
# =====================================================================
# Step 2: Architeuthis false-positive filter -> Kraken-style report -> Bracken
# One array task per sample; small memory. Submitted by submit_all.sh.
# =====================================================================
#SBATCH --job-name=smdb_filter_bracken
#SBATCH --output=%x_%A_%a.log
set -euo pipefail

source "${PIPE_DIR:-$(dirname "$0")}/config.sh"
conda_init

load_sample "$SLURM_ARRAY_TASK_ID"
KRAKEN_OUTPUT="$KRAKEN_OUT_DIR/${SAMPLE}_kraken.output"
FILTERED="$FILTER_OUT_DIR/${SAMPLE}_filtered.output"
FILTERED_REPORT="$FILTER_OUT_DIR/${SAMPLE}_filtered.kreport"

echo "[$(date)] Sample $SAMPLE on $(hostname)"

if [[ -f "$FILTER_OUT_DIR/${SAMPLE}.done" ]]; then
    echo "Already finished; skipping. (Delete $FILTER_OUT_DIR/${SAMPLE}.done to redo.)"
    exit 0
fi

# Kraken output may already be compressed from an earlier run
if [[ ! -s "$KRAKEN_OUTPUT" && -s "$KRAKEN_OUTPUT.gz" ]]; then
    echo "Decompressing earlier Kraken2 output"
    gunzip -k "$KRAKEN_OUTPUT.gz"
fi
[[ -s "$KRAKEN_OUTPUT" ]] || { echo "ERROR: missing $KRAKEN_OUTPUT (did step 1 finish?)"; exit 1; }

# --- 2a. Architeuthis: drop reads whose k-mers point to inconsistent lineages
if [[ ! -s "$FILTERED" && ! -s "$FILTERED.gz" ]]; then
    echo "[$(date)] Architeuthis filter"
    architeuthis mapping filter "$KRAKEN_OUTPUT" \
        --db "$DBDIR" --data-dir "$DB_TAXONOMY_DIR" \
        --out "$FILTERED.tmp"
    mv "$FILTERED.tmp" "$FILTERED"
fi
if [[ ! -s "$FILTERED" && -s "$FILTERED.gz" ]]; then
    gunzip -k "$FILTERED.gz"
fi

# --- 2b. Convert filtered per-read output to a Kraken report (what Bracken reads)
echo "[$(date)] kraken2-report"
"$KRAKEN2_REPORT_BIN" "$DB_TAXO" "$FILTERED" "$FILTERED_REPORT.tmp"
sed -i '/^$/d' "$FILTERED_REPORT.tmp"     # authors' workflow strips blank lines too
mv "$FILTERED_REPORT.tmp" "$FILTERED_REPORT"
[[ -s "$FILTERED_REPORT" ]] || { echo "ERROR: empty report $FILTERED_REPORT"; exit 1; }

# --- 2c. Bracken at each requested rank
for LVL in $BRACKEN_LEVELS; do
    OUT="$BRACKEN_OUT_DIR/${SAMPLE}_${LVL}.bracken"
    echo "[$(date)] Bracken level $LVL"
    bracken -d "$DBDIR" -i "$FILTERED_REPORT" -r "$READ_LEN" -l "$LVL" \
        -t "$BRACKEN_THRESHOLD" \
        -o "$OUT" \
        -w "$BRACKEN_OUT_DIR/${SAMPLE}_${LVL}.bracken.kreport"
done

# --- 2d. Shrink the big per-read files now that the reports exist
if [[ "$COMPRESS_INTERMEDIATES" == "true" ]]; then
    echo "[$(date)] Compressing per-read outputs"
    for f in "$KRAKEN_OUTPUT" "$FILTERED"; do
        if [[ -s "$f" ]]; then
            gzip -f "$f"
        fi
    done
fi

touch "$FILTER_OUT_DIR/${SAMPLE}.done"
echo "[$(date)] Done"
