# =====================================================================
# SoilMicrobeDB pipeline configuration
# Edit the lines marked EDIT, then run:  bash 00_setup.sh
# Every other script sources this file.
# =====================================================================

# --- EDIT: where your metagenomes live (the folder containing PH28metagenome/, SWa13metagenome_2/, ...)
READS_DIR="$HOME/frey_ph_metagenomes"

# --- EDIT: where results should go (needs space: Kraken2 per-read output is
#     roughly 1-2x the size of the uncompressed reads, per sample, until step 2 compresses it)
OUT_DIR="$HOME/frey_ph_metagenomes/frey_ph_smdb_results"

# --- Database location (the folder you ran `aws s3 sync` into)
DB_ROOT="$HOME/../shared/SoilMicrobeDB"

# --- Conda environment
CONDA_ENV="soilMicrobeDB"

# --- EDIT if needed: how your cluster makes `conda` available inside a batch job.
#     Common alternatives:  module load miniconda   /   module load anaconda3
conda_init() {
    #if ! command -v conda >/dev/null 2>&1; then
        module load anaconda/colsa 2>/dev/null || module load anaconda3 2>/dev/null || true
    #fi
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
}

# --- EDIT if needed: SLURM options added to every sbatch call (account, QOS ...)
#     e.g. SBATCH_COMMON="--account=frey_lab --qos=normal"
SBATCH_COMMON="--mail-user=eric.morrison@unh.edu --mail-type=END,FAIL"

# --- Kraken2 memory mode -------------------------------------------------------
# "ram"  : load the whole ~280 GB hash into memory. Fastest. Needs your 515 GB nodes.
# "mmap" : --memory-mapping. Runs on 257 GB nodes but can be MUCH slower, because the
#          hash is read from disk on demand. Copying to node-local disk helps (see below).
KRAKEN_MODE="ram"
KRAKEN_MEM="340G"          # for "ram" mode; use ~240G for "mmap" on 257 GB nodes
KRAKEN_CPUS=24
KRAKEN_TIME="12:00:00"
KRAKEN_PARTITION="shared"        # e.g. "bigmem" if 515 GB nodes are in their own partition
KRAKEN_CONSTRAINT=""       # e.g. "mem512" if your site tags big nodes with a feature
KRAKEN_MAX_PARALLEL=2      # how many Kraken2 jobs may run at once (each holds a big node)

# "mmap" mode only: copy the database to node-local disk first (needs ~285 GB free there)
MMAP_COPY_TO_LOCAL="false"

# --- Step 2 (Architeuthis + Bracken) resources --------------------------------
FILTER_MEM="32G"
FILTER_CPUS=4
FILTER_TIME="12:00:00"
FILTER_PARTITION=""

# --- Bracken -------------------------------------------------------------------
READ_LEN=150               # SoilMicrobeDB ships database150mers.kmer_distrib only
BRACKEN_THRESHOLD=10       # Bracken -t (min reads for a taxon to be re-estimated); 10 = default
BRACKEN_LEVELS="S G F O C P"   # species, genus, family, order, class, phylum

# --- Scratch for de-interleaved reads (node-local is best). $TMPDIR is set by many SLURM sites.
LOCAL_SCRATCH="${TMPDIR:-/tmp}"

# --- Housekeeping ----------------------------------------------------------------
# After step 2 succeeds, gzip the large per-read Kraken2 / Architeuthis outputs.
COMPRESS_INTERMEDIATES="true"

# =====================================================================
# Derived paths -- normally no need to edit below this line
# =====================================================================
PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_DIR="$PIPE_DIR/tools"
SAMPLES_TSV="$PIPE_DIR/samples.tsv"
KRAKEN2_REPORT_BIN="$TOOLS_DIR/kraken2-report"

# Accept either DB_ROOT itself or DB_ROOT/kraken2 as the folder holding hash.k2d
if [[ -f "$DB_ROOT/hash.k2d" ]]; then
    DBDIR="$DB_ROOT"
else
    DBDIR="$DB_ROOT/kraken2"
fi
DB_TAXONOMY_DIR="$DBDIR/taxonomy"
DB_TAXO="$DBDIR/taxo.k2d"

KRAKEN_OUT_DIR="$OUT_DIR/01_kraken"
FILTER_OUT_DIR="$OUT_DIR/02_architeuthis"
BRACKEN_OUT_DIR="$OUT_DIR/03_bracken"
SUMMARY_DIR="$OUT_DIR/04_summary"
LOG_DIR="$OUT_DIR/logs"

# Look up one line of samples.tsv by 1-based index: sets SAMPLE and FASTQ
load_sample() {
    local idx="$1"
    local line
    line="$(awk -F'\t' -v i="$idx" 'NR == i + 1 { print; exit }' "$SAMPLES_TSV")"
    if [[ -z "$line" ]]; then
        echo "ERROR: no sample at index $idx in $SAMPLES_TSV" >&2
        exit 1
    fi
    SAMPLE="$(cut -f1 <<<"$line")"
    FASTQ="$(cut -f2 <<<"$line")"
}
