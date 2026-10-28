# SoilMicrobeDB pipeline for the Frey pH metagenomes (SLURM)

Kraken2 -> Architeuthis -> Bracken, following the SoilMicrobeDB tutorial
(Werbin et al. 2025), adapted for JGI interleaved FASTQs and a SLURM cluster.

## Files

| File | What it does |
|---|---|
| `config.sh` | All paths and resources. **Edit `READS_DIR` and `OUT_DIR` first.** |
| `00_setup.sh` | One-time: installs missing tools into `soilMicrobeDB`, compiles `kraken2-report`, checks the DB, finds samples, checks reads. |
| `01_kraken2.slurm` | Array job, one task per sample: de-interleaves reads on node-local scratch, runs Kraken2. Big memory. |
| `02_filter_bracken.slurm` | Array job: Architeuthis filter -> `kraken2-report` -> Bracken at S/G/F/O/C/P. Small memory. |
| `03_summarize.slurm` | Combines samples into fungal summary tables. |
| `summarize_fungi.py` | The summary logic (standard-library Python). |
| `submit_all.sh` | Submits steps 1-3 with dependencies. |

## Run it

```bash
cd smdb_pipeline
nano config.sh                 # set READS_DIR and OUT_DIR (and partition/account if needed)
bash 00_setup.sh               # on a login node; fix anything marked [FAIL]
bash submit_all.sh 1           # test with one sample first
bash submit_all.sh             # then all 19 (finished samples are skipped)
```

Watch with `squeue -u $USER`; logs go to `$OUT_DIR/logs/`.

## Memory choice

`hash.k2d` is ~280 GB.

* **515 GB nodes (default, `KRAKEN_MODE="ram"`, `KRAKEN_MEM="340G"`)**: whole DB in RAM; fastest.
  If those nodes sit in a separate partition or have a feature tag, set `KRAKEN_PARTITION`
  or `KRAKEN_CONSTRAINT` (`sinfo -o "%P %m %f"` shows partitions, memory and features).
* **257 GB nodes (`KRAKEN_MODE="mmap"`, `KRAKEN_MEM="240G"`)**: uses `--memory-mapping`.
  Works, but may be many times slower because the hash is read from disk.
  Setting `MMAP_COPY_TO_LOCAL="true"` copies the DB to node-local disk first, which helps if
  the nodes have ~300 GB free local disk.
* 128 GB nodes are not a good fit for Kraken2 with this database; they're fine for step 2.

## Outputs (`$OUT_DIR`)

```
01_kraken/        <sample>_kraken.kreport, <sample>_kraken.output(.gz)
02_architeuthis/  <sample>_filtered.kreport, <sample>_filtered.output(.gz)
03_bracken/       <sample>_<RANK>.bracken  (+ .bracken.kreport)
04_summary/
  kingdom_summary.tsv                    fungal reads and % of total / classified / eukaryotic reads,
                                         after filtering, plus the unfiltered numbers for comparison
  fungi_<RANK>_reads.tsv                 fungal taxa x samples, Bracken-estimated reads
  fungi_<RANK>_relabund_within_fungi.tsv fungal community composition (each column sums to 1)
  fungi_<RANK>_share_of_all_taxa.tsv     fungal taxa as a share of everything Bracken estimated at that rank
  bracken_<RANK>_all_taxa_long.tsv       every taxon (bacteria, archaea, fungi...) in long format with lineage
```

For "how much of the soil community is fungal", use `kingdom_summary.tsv`.
For "which fungi, in what proportions", use `fungi_S_relabund_within_fungi.tsv` or the genus version.

## Notes

* Bracken uses `database150mers.kmer_distrib`; `00_setup.sh` warns if your reads aren't ~150 bp.
* Kraken2 per-read outputs are large; step 2 gzips them afterward (`COMPRESS_INTERMEDIATES`).
* Re-running is safe: finished samples are skipped at each step.
