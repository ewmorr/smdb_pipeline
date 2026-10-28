#!/bin/bash
# =====================================================================
# One-time setup + sanity checks. Run on a login node:   bash 00_setup.sh
#   - installs any missing tools into the conda env
#   - compiles kraken2-report (needed to turn Architeuthis output into a report)
#   - checks the database download is complete
#   - finds your samples and writes samples.tsv
#   - checks the reads are interleaved pairs and ~150 bp
# Safe to re-run.
# =====================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
conda_init

ok()   { echo "  [ok]   $*"; }
warn() { echo "  [WARN] $*"; }
fail() { echo "  [FAIL] $*"; FAILED=1; }
FAILED=0

echo "== 1. Software in conda env '$CONDA_ENV'"
need=()
command -v kraken2      >/dev/null || need+=("bioconda::kraken2")
command -v bracken      >/dev/null || need+=("bioconda::bracken")
command -v architeuthis >/dev/null || need+=("bioconda::architeuthis")
command -v taxonkit     >/dev/null || need+=("bioconda::taxonkit")
command -v reformat.sh  >/dev/null || need+=("bioconda::bbmap")
command -v python3      >/dev/null || need+=("conda-forge::python")
if (( ${#need[@]} )); then
    echo "  Installing: ${need[*]}"
    installer=conda; command -v mamba >/dev/null && installer=mamba
    "$installer" install -y -n "$CONDA_ENV" -c conda-forge -c bioconda "${need[@]}"
    conda activate "$CONDA_ENV"
fi
for t in kraken2 bracken architeuthis taxonkit reformat.sh python3; do
    if command -v "$t" >/dev/null; then ok "$t -> $(command -v "$t")"; else fail "$t not found"; fi
done

echo "== 2. kraken2-report"
mkdir -p "$TOOLS_DIR"
if [[ -x "$KRAKEN2_REPORT_BIN" ]]; then
    ok "already built: $KRAKEN2_REPORT_BIN"
else
    src="$TOOLS_DIR/kraken2_src"
    if [[ ! -d "$src" ]]; then
        # The SoilMicrobeDB repo carries the kraken2 source with kraken2-report.cpp added
        git clone --depth 1 https://github.com/zoey-rw/SoilMicrobeDB.git "$TOOLS_DIR/SoilMicrobeDB_repo"
        cp -r "$TOOLS_DIR/SoilMicrobeDB_repo/AWS_conda_environment_example/kraken2-master/src" "$src"
    fi
    CXX_BIN="${CXX:-}"
    [[ -z "$CXX_BIN" ]] && command -v g++ >/dev/null && CXX_BIN=g++
    if [[ -z "$CXX_BIN" ]]; then
        module load gcc 2>/dev/null || true
        command -v g++ >/dev/null && CXX_BIN=g++
    fi
    if [[ -z "$CXX_BIN" ]]; then
        echo "  No C++ compiler found; installing one into the env"
        conda install -y -n "$CONDA_ENV" -c conda-forge cxx-compiler
        conda activate "$CONDA_ENV"
        CXX_BIN="${CXX:-g++}"
    fi
    ( cd "$src" && "$CXX_BIN" -O3 -std=c++11 mmap_file.cc reports.cc taxonomy.cc kraken2-report.cpp \
          -o "$KRAKEN2_REPORT_BIN" )
    [[ -x "$KRAKEN2_REPORT_BIN" ]] && ok "built $KRAKEN2_REPORT_BIN" || fail "kraken2-report did not compile"
fi

echo "== 3. Database at $DBDIR"
for f in hash.k2d opts.k2d taxo.k2d "database${READ_LEN}mers.kmer_distrib" taxonomy/nodes.dmp taxonomy/names.dmp; do
    if [[ -s "$DBDIR/$f" ]]; then ok "$f ($(du -h "$DBDIR/$f" | cut -f1))"; else fail "missing $DBDIR/$f"; fi
done
if [[ -s "$DBDIR/hash.k2d" ]]; then
    hsize=$(stat -c %s "$DBDIR/hash.k2d")
    (( hsize < 250000000000 )) && warn "hash.k2d is only $((hsize/1000000000)) GB; expected ~280 GB. Is the download finished?"
fi

echo "== 4. Samples under $READS_DIR"
[[ -d "$READS_DIR" ]] || { fail "READS_DIR does not exist -- edit config.sh"; exit 1; }
printf "sample\tfastq\n" > "$SAMPLES_TSV"
n=0
while IFS= read -r fq; do
    # sample name = top-level folder, e.g. PH28metagenome_2
    rel="${fq#"$READS_DIR"/}"
    samp="${rel%%/*}"
    printf "%s\t%s\n" "$samp" "$fq" >> "$SAMPLES_TSV"
    n=$((n + 1))
done < <(find "$READS_DIR" -path "*/Filtered_Raw_Data/*.filter-METAGENOME.fastq.gz" -type f | sort -V)
dups=$(tail -n +2 "$SAMPLES_TSV" | cut -f1 | sort | uniq -d)
[[ -n "$dups" ]] && fail "more than one FASTQ for: $dups"
(( n > 0 )) && ok "$n samples written to $SAMPLES_TSV" || fail "no *.filter-METAGENOME.fastq.gz files found"

echo "== 5. Read format (first 20,000 reads of each file)"
while IFS=$'\t' read -r samp fq; do
    res=$( { zcat "$fq" 2>/dev/null || true; } | head -n 80000 | awk '
        NR % 4 == 1 { split($1, a, "/"); id = a[1]; sub(/^@/, "", id)
                      if ((NR - 1) % 8 == 0) { prev = id } else if (id != prev) { bad++ } }
        NR % 4 == 2 { n++; s += length($0); if (length($0) > mx) mx = length($0) }
        END { printf "%d %d %d", (n ? s / n : 0), mx, bad + 0 }')
    read -r mean mx bad <<<"$res"
    msg="$samp: mean ${mean} bp, max ${mx} bp"
    if (( bad > 0 )); then
        fail "$msg -- reads do NOT look interleaved ($bad unpaired names)"
    elif (( mx < READ_LEN - 10 || mx > READ_LEN + 10 )); then
        warn "$msg -- max length differs from READ_LEN=$READ_LEN; Bracken may need a matching kmer_distrib"
    else
        ok "$msg, interleaved pairs"
    fi
done < <(tail -n +2 "$SAMPLES_TSV")

echo "== 6. Output folders under $OUT_DIR"
mkdir -p "$KRAKEN_OUT_DIR" "$FILTER_OUT_DIR" "$BRACKEN_OUT_DIR" "$SUMMARY_DIR" "$LOG_DIR"
ok "created"
echo "  Free space there: $(df -h "$OUT_DIR" | awk 'NR==2 {print $4}')"

echo
if (( FAILED )); then
    echo "Setup found problems (see [FAIL] above). Fix them, then re-run."
    exit 1
fi
echo "Setup complete. Next:  bash submit_all.sh"
