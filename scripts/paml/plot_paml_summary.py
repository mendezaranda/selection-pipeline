#!/usr/bin/env python3
"""
2-panel summary figure for a codeml branch-model scan, for any phenotype:
  A) p-value histogram (diagnostic for genome-wide scan behavior)
  B) "Manhattan-style" scan across orthogroups, significance threshold line
     (BH-FDR or Holm-Bonferroni, your choice), top hits labeled by gene name
     (falls back to OG id if unannotated)

Takes --phenotype instead of an output filename -- the phenotype name plus
the chosen correction drive the output filename automatically, e.g.
--phenotype acid --correction BH -> acid_summary_BH.png.

Usage:
    python3 plot_paml_summary.py <annotated_significance.tsv> --phenotype acid \
        [--correction BH|holm] [--top-n 8] [--title "..."] [--outdir .]

<annotated_significance.tsv>: output of annotate_with_gene_names.py run on a
    build_significance_table.py output -- needs OG, gene_symbol, p_raw, and
    q_BH and/or p_holm (whichever --correction you pick).
"""
import sys
import os
import csv
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# correction -> (column name in the input tsv, human-readable label, output filename tag)
CORRECTIONS = {
    "bh":   {"col": "q_BH",   "label": "BH-FDR q < 0.05",          "thresh_label": "BH-FDR 5%",          "tag": "BH"},
    "holm": {"col": "p_holm", "label": "Holm-Bonferroni p < 0.05", "thresh_label": "Holm-Bonferroni 5%", "tag": "Holm"},
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_tsv")
    ap.add_argument("--phenotype", required=True,
                     help="Phenotype name, e.g. acid -- drives the output filename "
                          "(<phenotype>_summary_<correction>.png) and the default title.")
    ap.add_argument("--correction", default="BH", choices=["BH", "bh", "holm", "Holm", "HOLM"],
                     help="Which multiple-testing correction defines 'significant' for the "
                          "threshold line, point coloring, and the printed top-10 ranking. "
                          "Default: BH (Benjamini-Hochberg FDR).")
    ap.add_argument("--title", default=None,
                     help="Panel B title. Default: '<phenotype> foreground'.")
    ap.add_argument("--top-n", type=int, default=8)
    ap.add_argument("--outdir", default=".", help="Directory to write the PNG into. Default: current dir.")
    args = ap.parse_args()

    corr = CORRECTIONS[args.correction.lower()]
    title = args.title if args.title is not None else f"{args.phenotype} foreground"
    output_png = os.path.join(args.outdir, f"{args.phenotype}_summary_{corr['tag']}.png")

    with open(args.input_tsv) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))

    if rows and corr["col"] not in rows[0]:
        sys.exit(f"'{args.input_tsv}' has no '{corr['col']}' column (needed for "
                 f"--correction {args.correction}). Found columns: {list(rows[0].keys())}")

    n = len(rows)
    p = np.array([float(r['p_raw']) for r in rows])
    qc = np.array([float(r[corr["col"]]) for r in rows])
    genes = np.array([r.get('gene_symbol', 'NA') for r in rows])
    ogs = np.array([r['OG'] for r in rows])
    neglog10p = -np.log10(np.clip(p, 1e-300, 1))
    sig = qc < 0.05

    COL_BG = "#8A94A6"     # muted slate, non-significant
    COL_SIG = "#D6604D"    # warm red-orange, significant
    COL_LINE = "#3D4451"

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.0), dpi=200)

    # --- Panel A: p-value histogram ---
    ax = axes[0]
    ax.hist(p, bins=40, color=COL_BG, edgecolor="white", linewidth=0.4)
    ax.set_xlabel("raw p-value (branch-model LRT, df=1)")
    ax.set_ylabel("number of orthogroups")
    ax.set_title("A. p-value distribution", loc="left", fontsize=11, fontweight="bold")
    ax.spines[['top', 'right']].set_visible(False)
    ax.text(0.98, 0.95, f"n = {n:,} orthogroups tested",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color="#555")

    # --- Panel B: Manhattan-style scan, top hits labeled ---
    ax = axes[1]
    order = np.argsort(ogs)
    x = np.arange(n)
    sig_o = sig[order]
    ax.scatter(x[~sig_o], neglog10p[order][~sig_o], s=6, color=COL_BG, alpha=0.6,
               linewidths=0, label="not significant")
    ax.scatter(x[sig_o], neglog10p[order][sig_o], s=10, color=COL_SIG, alpha=0.9,
               linewidths=0, label=corr["label"])

    sig_ps = p[sig]
    if len(sig_ps):
        thresh = sig_ps.max()
        ax.axhline(-np.log10(thresh), color=COL_LINE, linestyle="--", linewidth=1, alpha=0.8)
        ax.text(x[-1], -np.log10(thresh), f"  {corr['thresh_label']}", va="bottom", ha="right",
                fontsize=8, color=COL_LINE)

    pos_in_order = {orig_i: xi for xi, orig_i in enumerate(order)}
    top_idx_global = np.argsort(p)[:args.top_n]
    for gi in top_idx_global:
        xi = pos_in_order[gi]
        yi = neglog10p[gi]
        label = genes[gi] if genes[gi] not in ("NA", "") else ogs[gi]
        ax.annotate(label, (xi, yi), textcoords="offset points", xytext=(3, 4),
                    fontsize=7.5, color=COL_LINE, fontstyle="italic")

    ax.set_xlabel("orthogroup (ordered by OG ID)")
    ax.set_ylabel("-log10(p)")
    ax.set_title(f"B. {title}", loc="left", fontsize=11, fontweight="bold")
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    plt.tight_layout()
    if args.outdir != ".":
        os.makedirs(args.outdir, exist_ok=True)
    plt.savefig(output_png, facecolor="white")
    print(f"Saved {output_png}  (correction: {corr['tag']}, {int(sig.sum())}/{n} significant)")

    top10 = sorted(rows, key=lambda r: float(r[corr["col"]]))[:10]
    print(f"\nTop 10 by {corr['col']}:")
    for r in top10:
        print(f"  {r['OG']:<12}{r.get('gene_symbol','NA'):<14}{corr['col']}={float(r[corr['col']]):.4g}"
              f"  p_raw={float(r['p_raw']):.4g}")


if __name__ == "__main__":
    main()
