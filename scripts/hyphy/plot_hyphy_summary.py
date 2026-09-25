#!/usr/bin/env python3
"""
2-panel summary figure for a HyPhy RELAX or aBSREL genome-wide scan, for any
phenotype -- same layout/logic as plot_paml_summary.py:
  A) p-value histogram (diagnostic for genome-wide scan behavior)
  B) "Manhattan-style" scan across orthogroups, significance threshold line
     (BH-FDR or Holm-Bonferroni, your choice), top hits labeled by gene name
     (falls back to OG id if unannotated)

Takes --phenotype and --test instead of an output filename -- these plus the
chosen correction drive the output filename automatically, e.g.
--phenotype acid --test absrel --correction BH -> acid_hyphy_absrel_summary_BH.png.

RELAX and aBSREL are different hypothesis tests (relaxed/intensified
selection strength vs. episodic diversifying selection on the Test
branches) with their own p-value/correction columns in the significance
table (from build_hyphy_significance_table.py) -- --test picks which one to
plot; run this script twice (once per --test) to get both figures.

Usage:
    python3 plot_hyphy_summary.py <annotated_significance.tsv> --phenotype acid \
        --test relax|absrel [--correction BH|holm] [--top-n 8] [--title "..."] [--outdir .]

<annotated_significance.tsv>: output of annotate_with_gene_names.py run on a
    build_hyphy_significance_table.py output -- needs OG, gene_symbol, plus
    RELAX_p/RELAX_q_BH/RELAX_p_holm (--test relax) or
    ABSREL_min_p/ABSREL_q_BH/ABSREL_p_holm (--test absrel).
"""
import sys
import os
import csv
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# test -> (raw p-value column, human label, x-axis label for panel A)
TESTS = {
    "relax": {
        "p_col": "RELAX_p",
        "label": "RELAX",
        "xlabel_a": "raw p-value (RELAX LR test)",
        "q_bh_col": "RELAX_q_BH",
        "p_holm_col": "RELAX_p_holm",
    },
    "absrel": {
        "p_col": "ABSREL_min_p",
        "label": "aBSREL",
        "xlabel_a": "min corrected p-value across Test branches (aBSREL)",
        "q_bh_col": "ABSREL_q_BH",
        "p_holm_col": "ABSREL_p_holm",
    },
}

# correction -> (label, threshold label, output filename tag)
CORRECTIONS = {
    "bh":   {"label": "BH-FDR q < 0.05",          "thresh_label": "BH-FDR 5%",          "tag": "BH"},
    "holm": {"label": "Holm-Bonferroni p < 0.05", "thresh_label": "Holm-Bonferroni 5%", "tag": "Holm"},
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_tsv")
    ap.add_argument("--phenotype", required=True,
                     help="Phenotype name, e.g. acid -- drives the output filename "
                          "(<phenotype>_hyphy_<test>_summary_<correction>.png) and the default title.")
    ap.add_argument("--test", required=True, choices=["relax", "absrel", "RELAX", "ABSREL"],
                     help="Which HyPhy test to plot: relax or absrel.")
    ap.add_argument("--correction", default="BH", choices=["BH", "bh", "holm", "Holm", "HOLM"],
                     help="Which multiple-testing correction defines 'significant' for the "
                          "threshold line, point coloring, and the printed top-10 ranking. "
                          "Default: BH (Benjamini-Hochberg FDR).")
    ap.add_argument("--title", default=None,
                     help="Panel B title. Default: '<phenotype> foreground (<test>)'.")
    ap.add_argument("--top-n", type=int, default=8)
    ap.add_argument("--outdir", default=".", help="Directory to write the PNG into. Default: current dir.")
    args = ap.parse_args()

    test = TESTS[args.test.lower()]
    corr = CORRECTIONS[args.correction.lower()]
    q_col = test["q_bh_col"] if corr["tag"] == "BH" else test["p_holm_col"]
    title = args.title if args.title is not None else f"{args.phenotype} foreground ({test['label']})"
    output_png = os.path.join(args.outdir, f"{args.phenotype}_hyphy_{args.test.lower()}_summary_{corr['tag']}.png")

    with open(args.input_tsv) as f:
        all_rows = list(csv.DictReader(f, delimiter='\t'))

    if all_rows and (test["p_col"] not in all_rows[0] or q_col not in all_rows[0]):
        sys.exit(f"'{args.input_tsv}' is missing '{test['p_col']}' or '{q_col}' "
                 f"(needed for --test {args.test} --correction {args.correction}). "
                 f"Found columns: {list(all_rows[0].keys())}")

    # Only rows where this specific test actually ran (build_hyphy_significance_table.py
    # writes 'NA' for OGs where this test was skipped -- e.g. tree-pruning left no Test
    # taxon for that OG, or you haven't run this test yet with ONLY=).
    rows = [r for r in all_rows if r[test["p_col"]] not in ("NA", "")]
    n_skipped = len(all_rows) - len(rows)
    if not rows:
        sys.exit(f"No rows with a valid {test['p_col']} in '{args.input_tsv}' -- "
                 f"has {test['label']} actually been run for this phenotype yet?")

    n = len(rows)
    p = np.array([float(r[test["p_col"]]) for r in rows])
    qc = np.array([float(r[q_col]) for r in rows])
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
    ax.set_xlabel(test["xlabel_a"])
    ax.set_ylabel("number of orthogroups")
    ax.set_title("A. p-value distribution", loc="left", fontsize=11, fontweight="bold")
    ax.spines[['top', 'right']].set_visible(False)
    note = f"n = {n:,} orthogroups tested"
    if n_skipped:
        note += f"\n({n_skipped:,} skipped: no {test['label']} result for this OG)"
    ax.text(0.98, 0.95, note, transform=ax.transAxes, ha="right", va="top", fontsize=8, color="#555")

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
    print(f"Saved {output_png}  ({test['label']}, correction: {corr['tag']}, "
          f"{int(sig.sum())}/{n} significant, {n_skipped} OGs had no {test['label']} result)")

    top10 = sorted(rows, key=lambda r: float(r[test["p_col"]]))[:10]
    print(f"\nTop 10 by {test['p_col']}:")
    for r in top10:
        print(f"  {r['OG']:<12}{r.get('gene_symbol','NA'):<14}{q_col}={float(r[q_col]):.4g}"
              f"  {test['p_col']}={float(r[test['p_col']]):.4g}")


if __name__ == "__main__":
    main()
