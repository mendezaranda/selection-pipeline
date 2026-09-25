#!/usr/bin/env python3
"""
Builds a ranked significance table (BH-FDR + Holm-Bonferroni corrections)
from parse_hyphy_results.py's combined RELAX+aBSREL TSV. Same corrections,
same conventions as build_significance_table.py for codeml -- but RELAX and
aBSREL are two DIFFERENT hypothesis tests (relaxed/intensified selection
strength vs. episodic diversifying selection), each with its own p-value, so
each gets corrected independently rather than picking one "the" p-value.
An OG missing one test's result (e.g. you ran ONLY=absrel and haven't run
RELAX yet, or a test was skipped for that OG by the tree-pruning check in
run_hyphy_relax_absrel.sh) is simply excluded from THAT test's correction
and ranking, not dropped from the table entirely.

Input columns expected (as produced by parse_hyphy_results.py):
    OG  RELAX_K  RELAX_LR  RELAX_p  ABSREL_min_p  ABSREL_min_p_branch
    ABSREL_n_test_branches_p05  ABSREL_n_test_branches_total

Usage:
    python3 build_hyphy_significance_table.py <hyphy_results.tsv> <output_significance.tsv>
"""
import sys
import csv


def bh_fdr(rows, p_key, q_key):
    """In-place: adds q_key (BH-FDR q-value) to every dict in rows that has
    a numeric p_key. rows must already be sorted ascending by p_key."""
    n = len(rows)
    prev_q = 1.0
    for i, r in enumerate(reversed(rows)):
        rank = n - i
        q = r[p_key] * n / rank
        q = min(q, prev_q)
        prev_q = q
        r[q_key] = q


def holm(rows, p_key, out_key):
    """In-place: adds out_key (Holm-Bonferroni corrected p) to every dict in
    rows. rows must already be sorted ascending by p_key."""
    n = len(rows)
    prev_holm = 0.0
    for i, r in enumerate(rows):
        rank = i + 1
        h = r[p_key] * (n - rank + 1)
        h = max(h, prev_holm)
        h = min(h, 1.0)
        prev_holm = h
        r[out_key] = h


def correct_one_test(all_rows, p_key, q_out, holm_out, label):
    """Selects rows with a valid numeric p_key, BH+Holm-corrects them, and
    writes the results back onto the SAME dicts in all_rows (rows without a
    valid p_key get 'NA' for both correction columns). Returns the count of
    testable rows and how many pass BH q<0.05 / Holm p<0.05, for the summary
    printout."""
    testable = []
    for r in all_rows:
        raw = r.get(p_key)
        if raw in (None, "NA", ""):
            r[q_out] = "NA"
            r[holm_out] = "NA"
            continue
        try:
            r[p_key + "_float"] = float(raw)
        except ValueError:
            r[q_out] = "NA"
            r[holm_out] = "NA"
            continue
        testable.append(r)

    testable.sort(key=lambda r: r[p_key + "_float"])
    bh_fdr(testable, p_key + "_float", q_out)
    holm(testable, p_key + "_float", holm_out)
    for r in testable:
        del r[p_key + "_float"]

    n = len(testable)
    sig_raw = sum(1 for r in testable if float(r[p_key]) < 0.05)
    sig_bh = sum(1 for r in testable if r[q_out] != "NA" and r[q_out] < 0.05)
    sig_holm = sum(1 for r in testable if r[holm_out] != "NA" and r[holm_out] < 0.05)
    print(f"{label}: {n} testable OGs -- raw p<0.05: {sig_raw} "
          f"({100*sig_raw/n:.1f}%)  BH q<0.05: {sig_bh}  Holm p<0.05: {sig_holm}"
          if n else f"{label}: 0 testable OGs")
    return n


def main():
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <hyphy_results.tsv> <output_significance.tsv>")
    in_path, out_path = sys.argv[1:3]

    with open(in_path) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))

    print(f"Input rows (OGs with at least a RELAX or aBSREL result): {len(rows)}")
    correct_one_test(rows, "RELAX_p", "RELAX_q_BH", "RELAX_p_holm", "RELAX")
    correct_one_test(rows, "ABSREL_min_p", "ABSREL_q_BH", "ABSREL_p_holm", "aBSREL")

    fieldnames = [
        "OG", "RELAX_K", "RELAX_LR", "RELAX_p", "RELAX_q_BH", "RELAX_p_holm",
        "ABSREL_min_p", "ABSREL_min_p_branch", "ABSREL_q_BH", "ABSREL_p_holm",
        "ABSREL_n_test_branches_p05", "ABSREL_n_test_branches_total",
    ]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda r: r["OG"]):
            for key in fieldnames:
                r.setdefault(key, "NA")
            w.writerow(r)

    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
