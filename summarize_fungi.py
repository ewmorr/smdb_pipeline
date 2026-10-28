#!/usr/bin/env python3
"""
Summarize SoilMicrobeDB results across samples, with a focus on fungi.

Reads (from the pipeline's output folders):
  01_kraken/<sample>_kraken.kreport           raw Kraken2 report (full taxonomy, zero counts included)
  02_architeuthis/<sample>_filtered.kreport   report after Architeuthis filtering
  03_bracken/<sample>_<LEVEL>.bracken         Bracken abundance tables

Writes to 04_summary/:
  kingdom_summary.tsv                 per sample: reads in Bacteria / Archaea / Eukaryota / Fungi,
                                      and fungal share of all and of classified reads
                                      (before and after Architeuthis filtering)
  bracken_<L>_all_taxa_long.tsv       every Bracken taxon, every sample, with lineage + fungal flag
  fungi_<L>_reads.tsv                 fungal taxa x samples, Bracken estimated reads
  fungi_<L>_relabund_within_fungi.tsv fungal taxa x samples, share of fungal reads (columns sum to 1)
  fungi_<L>_share_of_all_taxa.tsv       fungal taxa x samples, share of ALL taxa (bacteria too) that Bracken estimated at that rank

Uses only the Python standard library.
"""
import argparse
import csv
import glob
import os
import sys
from collections import defaultdict

FUNGI_TAXID = "4751"   # NCBI taxonomy ID for kingdom Fungi
RANK_ORDER = ["D", "K", "P", "C", "O", "F", "G", "S"]


def read_kreport(path):
    """Yield (clade_reads, rank, taxid, depth, name) for each line of a Kraken report.

    Handles both the standard 6-column report and the 8-column
    --report-minimizer-data variant: rank, taxid and name are always the last three columns.
    """
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 6:
                continue
            raw_name = cols[-1]
            name = raw_name.lstrip(" ")
            depth = (len(raw_name) - len(name)) // 2
            yield int(cols[1]), cols[-3].strip(), cols[-2].strip(), depth, name


def build_lineages(kreport_paths):
    """taxid -> dict(name, rank, in_fungi, domain, lineage)."""
    info = {}
    for path in kreport_paths:
        stack = []  # (depth, taxid, rank, name)
        for _, rank, taxid, depth, name in read_kreport(path):
            if taxid in info:
                # still need the stack to stay correct for children
                while stack and stack[-1][0] >= depth:
                    stack.pop()
                stack.append((depth, taxid, rank, name))
                continue
            while stack and stack[-1][0] >= depth:
                stack.pop()
            stack.append((depth, taxid, rank, name))
            chain = stack
            in_fungi = any(t == FUNGI_TAXID or (r == "K" and n == "Fungi") for _, t, r, n in chain)
            domain = next((n for _, _, r, n in chain if r == "D"), "")
            lineage = ";".join(f"{r.lower()}__{n}" for _, _, r, n in chain if r in RANK_ORDER)
            info[taxid] = {"name": name, "rank": rank, "in_fungi": in_fungi,
                           "domain": domain, "lineage": lineage}
    return info


def clade_counts(path):
    """Return dict with unclassified, total, and clade reads for key groups."""
    out = defaultdict(int)
    for clade, rank, taxid, depth, name in read_kreport(path):
        if rank == "U":
            out["unclassified"] += clade
        elif rank == "R" and depth == 0:
            out["classified"] += clade
        elif rank == "D" and name in ("Bacteria", "Archaea", "Eukaryota", "Viruses"):
            out[name] += clade
        if taxid == FUNGI_TAXID or (rank == "K" and name == "Fungi"):
            out["Fungi"] = clade
    out["total"] = out["unclassified"] + out["classified"]
    return out


def pct(a, b):
    return f"{100.0 * a / b:.4f}" if b else "NA"


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True, help="pipeline OUT_DIR")
    ap.add_argument("--samples", required=True, help="samples.tsv")
    ap.add_argument("--levels", default="S G F O C P")
    args = ap.parse_args()

    kdir = os.path.join(args.out_dir, "01_kraken")
    fdir = os.path.join(args.out_dir, "02_architeuthis")
    bdir = os.path.join(args.out_dir, "03_bracken")
    sdir = os.path.join(args.out_dir, "04_summary")
    os.makedirs(sdir, exist_ok=True)

    with open(args.samples) as fh:
        samples = [r["sample"] for r in csv.DictReader(fh, delimiter="\t")]

    raw_reports = {s: os.path.join(kdir, f"{s}_kraken.kreport") for s in samples}
    filt_reports = {s: os.path.join(fdir, f"{s}_filtered.kreport") for s in samples}
    missing = [s for s in samples if not os.path.exists(filt_reports[s])]
    if missing:
        print(f"WARNING: no filtered report yet for {len(missing)} sample(s): {', '.join(missing)}",
              file=sys.stderr)
    samples = [s for s in samples if s not in missing]
    if not samples:
        sys.exit("No finished samples to summarize.")

    # Raw reports carry the full taxonomy (--report-zero-counts), so they give lineages for every taxid
    lineage_sources = [p for p in raw_reports.values() if os.path.exists(p)] + list(filt_reports.values())
    tax = build_lineages(lineage_sources)

    # ---- kingdom-level summary
    rows = []
    for s in samples:
        f = clade_counts(filt_reports[s])
        r = clade_counts(raw_reports[s]) if os.path.exists(raw_reports[s]) else defaultdict(int)
        rows.append([
            s, f["total"], f["classified"], f["Bacteria"], f["Archaea"], f["Eukaryota"], f["Fungi"],
            pct(f["Fungi"], f["total"]), pct(f["Fungi"], f["classified"]), pct(f["Fungi"], f["Eukaryota"]),
            r["classified"], r["Fungi"], pct(r["Fungi"], r["total"]),
        ])
    write_tsv(os.path.join(sdir, "kingdom_summary.tsv"),
              ["sample", "total_reads", "classified_reads", "bacteria_reads", "archaea_reads",
               "eukaryota_reads", "fungi_reads", "fungi_pct_of_total", "fungi_pct_of_classified",
               "fungi_pct_of_eukaryota", "unfiltered_classified_reads", "unfiltered_fungi_reads",
               "unfiltered_fungi_pct_of_total"], rows)

    # ---- Bracken tables per rank
    for lvl in args.levels.split():
        long_rows = []
        fungi_reads = defaultdict(dict)    # taxid -> sample -> reads
        fungi_frac_all = defaultdict(dict)
        found_any = False
        for s in samples:
            path = os.path.join(bdir, f"{s}_{lvl}.bracken")
            if not os.path.exists(path):
                print(f"WARNING: missing {path}", file=sys.stderr)
                continue
            found_any = True
            with open(path) as fh:
                for rec in csv.DictReader(fh, delimiter="\t"):
                    tid = rec["taxonomy_id"]
                    t = tax.get(tid, {"in_fungi": False, "domain": "", "lineage": ""})
                    reads = int(rec["new_est_reads"])
                    frac = float(rec["fraction_total_reads"])
                    long_rows.append([s, tid, rec["name"], t["domain"], int(t["in_fungi"]),
                                      t["lineage"], reads, f"{frac:.6g}"])
                    if t["in_fungi"]:
                        fungi_reads[tid][s] = reads
                        fungi_frac_all[tid][s] = frac
        if not found_any:
            continue
        write_tsv(os.path.join(sdir, f"bracken_{lvl}_all_taxa_long.tsv"),
                  ["sample", "taxid", "name", "domain", "is_fungus", "lineage",
                   "bracken_est_reads", "fraction_of_bracken_reads_at_rank"], long_rows)

        # wide fungal matrices, taxa sorted by total reads
        taxa = sorted(fungi_reads, key=lambda t: -sum(fungi_reads[t].values()))
        fungal_totals = {s: sum(fungi_reads[t].get(s, 0) for t in taxa) for s in samples}
        head = ["taxid", "name", "lineage"] + samples

        def wide(getter):
            return [[t, tax.get(t, {}).get("name", ""), tax.get(t, {}).get("lineage", "")]
                    + [getter(t, s) for s in samples] for t in taxa]

        write_tsv(os.path.join(sdir, f"fungi_{lvl}_reads.tsv"), head,
                  wide(lambda t, s: fungi_reads[t].get(s, 0)))
        write_tsv(os.path.join(sdir, f"fungi_{lvl}_relabund_within_fungi.tsv"), head,
                  wide(lambda t, s: f"{fungi_reads[t].get(s, 0) / fungal_totals[s]:.6g}"
                       if fungal_totals[s] else "0"))
        write_tsv(os.path.join(sdir, f"fungi_{lvl}_share_of_all_taxa.tsv"), head,
                  wide(lambda t, s: f"{fungi_frac_all[t].get(s, 0.0):.6g}"))
        print(f"Level {lvl}: {len(taxa)} fungal taxa across {len(samples)} samples")

    print(f"Summary written to {sdir}")


if __name__ == "__main__":
    main()
