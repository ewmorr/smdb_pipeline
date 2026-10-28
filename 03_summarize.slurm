#!/bin/bash
# =====================================================================
# Step 3: combine all samples into summary tables (fungal focus).
# Runs once after every step-2 task finishes. Also fine to run by hand:
#   bash 03_summarize.slurm
# =====================================================================
#SBATCH --job-name=smdb_summary
#SBATCH --output=%x_%j.log
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
set -euo pipefail

source "${PIPE_DIR:-$(dirname "$0")}/config.sh"
conda_init

python3 "$PIPE_DIR/summarize_fungi.py" \
    --out-dir "$OUT_DIR" \
    --samples "$SAMPLES_TSV" \
    --levels "$BRACKEN_LEVELS"
