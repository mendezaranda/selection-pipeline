#!/usr/bin/env python3
"""
Builds a ranked significance table (BH-FDR + Holm-Bonferroni corrections)
from a codeml branch-model summary TSV. Works for any phenotype -- just
point it at that phenotype's merged/deduped summary file.

Only rows with status == 'OK' and a valid p-value are tested/corrected;
everything else (SKIP_*, CODEML_FAIL) is dropped.

Input columns expected (as produced by og_to_codon_paml_codeml.sh /
codeml_from_existing_alignment.sh and merged the same way as before):
    OG  status  ntaxa  nsites  lnL_null  lnL_alt  LRT  p_df1  note

Usage:
    python3 build_significance_table.py <codeml_summary.deduped.tsv> <output_significance.tsv>

If your summary file is still the raw per-task/append-only version (not yet
deduped), dedupe it first the same way as before -- keep the LAST row per OG:
    awk -F'\t' 'NR==1{h=$0;next}{row[$1]=$0}END{print h; for (og in row) print row[og]}' \
        codeml_branch_models_summary.tsv | sort -k1,1 > codeml_branch_models_summary.deduped.tsv
"""
import sys
import csv


def main():
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <codeml_summary.deduped.tsv> <output_significance.tsv>")
    in_path, out_path = sys.argv[1:3]

    rows = []
    with open(in_path) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            if row['status'] != 'OK':
                continue
            p = row['p_df1']
            if p in ('NA', ''):
                continue
            rows.append({
                'OG': row['OG'],
                'ntaxa': int(row['ntaxa']),
                'nsites': int(row['nsites']),
                'lnL_null': float(row['lnL_null']),
                'lnL_alt': float(row['lnL_alt']),
                'LRT': float(row['LRT']),
                'p': float(p),
            })

    print(f"Testable OGs (status=OK, valid p-value): {len(rows)}")

    # Benjamini-Hochberg FDR
    rows_sorted = sorted(rows, key=lambda x: x['p'])
    n = len(rows_sorted)
    prev_q = 1.0
    for i, r in enumerate(reversed(rows_sorted)):
        rank = n - i
        q = r['p'] * n / rank
        q = min(q, prev_q)
        prev_q = q
        r['q_BH'] = q

    # Holm-Bonferroni
    prev_holm = 0.0
    for i, r in enumerate(rows_sorted):
        rank = i + 1
        holm = r['p'] * (n - rank + 1)
        holm = max(holm, prev_holm)
        holm = min(holm, 1.0)
        prev_holm = holm
        r['p_holm'] = holm

    sig_raw_05 = sum(1 for r in rows_sorted if r['p'] < 0.05)
    sig_bh_05 = sum(1 for r in rows_sorted if r['q_BH'] < 0.05)
    sig_bh_10 = sum(1 for r in rows_sorted if r['q_BH'] < 0.10)
    sig_holm_05 = sum(1 for r in rows_sorted if r['p_holm'] < 0.05)

    print(f"Raw p<0.05 (uncorrected):  {sig_raw_05} ({100*sig_raw_05/n:.1f}%)")
    print(f"BH-FDR q<0.05:             {sig_bh_05}")
    print(f"BH-FDR q<0.10:             {sig_bh_10}")
    print(f"Holm-Bonferroni p<0.05:    {sig_holm_05}")

    with open(out_path, "w") as o:
        o.write("OG\tntaxa\tnsites\tLRT\tp_raw\tq_BH\tp_holm\n")
        for r in rows_sorted:
            o.write(f"{r['OG']}\t{r['ntaxa']}\t{r['nsites']}\t{r['LRT']:.4f}\t"
                    f"{r['p']:.6g}\t{r['q_BH']:.6g}\t{r['p_holm']:.6g}\n")

    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
