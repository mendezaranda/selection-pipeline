#!/usr/bin/env python3
"""
Parses HyPhy RELAX + aBSREL JSON output (from run_hyphy_relax_absrel.sh)
into one combined TSV, one row per OG.

RELAX: reads the top-level "test results" block --
  LR, p-value, and the relaxation/intensification parameter K
  (K<1 = relaxed selection on Test branches, K>1 = intensified, relative to
  Reference).

aBSREL: reads "branch attributes" for just the Test-labeled branches (the
  ones run_hyphy_relax_absrel.sh restricted testing to via --branches Test)
  -- reports the minimum corrected p-value across those branches (i.e.
  "was at least one of the phenotype branches under episodic diversifying
  selection"), plus which branch(es) hit p<0.05.

This targets the HyPhy 2.5.x JSON schema. If your HyPhy version's JSON looks
different, this script will print which OGs it couldn't parse and why --
send me one raw <OG>.RELAX.json or <OG>.ABSREL.json and I'll fix the key
paths rather than guessing blind.

Usage:
    python3 parse_hyphy_results.py <hyphy_outdir> <absrel_fg_tree.nwk> <output.tsv>

<hyphy_outdir>: the OUTDIR passed to run_hyphy_relax_absrel.sh (needs
    relax/<OG>.RELAX.json and absrel/<OG>.ABSREL.json inside it).
<absrel_fg_tree.nwk>: the tree used for the aBSREL run (from
    make_hyphy_trees.py), to know which branch names were {Test}-labeled.
"""
import sys
import os
import re
import json
import glob


def get_test_branches(tree_path):
    text = open(tree_path).read()
    return set(re.findall(r"([A-Za-z][A-Za-z0-9_]*)\{Test\}", text))


def parse_relax(path):
    try:
        d = json.load(open(path))
    except Exception as e:
        return None, f"json_parse_error: {e}"
    tr = d.get("test results")
    if tr is None:
        return None, "no 'test results' key"
    try:
        return {
            "LR": tr["LR"],
            "p_value": tr["p-value"],
            "K": tr.get("relaxation or intensification parameter"),
        }, None
    except KeyError as e:
        return None, f"missing key {e}"


def parse_absrel(path, test_branches):
    try:
        d = json.load(open(path))
    except Exception as e:
        return None, f"json_parse_error: {e}"
    ba = d.get("branch attributes")
    if ba is None:
        return None, "no 'branch attributes' key"
    # partitions are keyed "0", "1", ... -- merge across all of them
    branch_p = {}
    for _part, branches in ba.items():
        if not isinstance(branches, dict):
            continue
        for bname, attrs in branches.items():
            if bname not in test_branches:
                continue
            if not isinstance(attrs, dict):
                continue
            p = attrs.get("Corrected P-value", attrs.get("Uncorrected P-value"))
            if p is not None:
                branch_p[bname] = p
    if not branch_p:
        return None, f"no p-values found for Test branches {sorted(test_branches)}"
    min_branch = min(branch_p, key=branch_p.get)
    return {
        "min_p": branch_p[min_branch],
        "min_p_branch": min_branch,
        "n_test_branches_p05": sum(1 for p in branch_p.values() if p < 0.05),
        "n_test_branches_total": len(branch_p),
        "per_branch": branch_p,
    }, None


def main():
    if len(sys.argv) != 4:
        sys.exit(f"usage: {sys.argv[0]} <hyphy_outdir> <absrel_fg_tree.nwk> <output.tsv>")
    outdir, absrel_tree, out_path = sys.argv[1:4]

    test_branches = get_test_branches(absrel_tree)
    if not test_branches:
        sys.exit(f"No {{Test}}-labeled branches found in {absrel_tree}")
    print(f"Test (foreground) branches: {sorted(test_branches)}")

    relax_files = sorted(glob.glob(os.path.join(outdir, "relax", "*.RELAX.json")))
    absrel_files = sorted(glob.glob(os.path.join(outdir, "absrel", "*.ABSREL.json")))
    print(f"Found {len(relax_files)} RELAX JSONs, {len(absrel_files)} aBSREL JSONs")

    ogs = set()
    relax_by_og = {}
    absrel_by_og = {}
    relax_errors = []
    absrel_errors = []

    for f in relax_files:
        og = os.path.basename(f).split(".")[0]
        ogs.add(og)
        res, err = parse_relax(f)
        if res:
            relax_by_og[og] = res
        else:
            relax_errors.append((og, err))

    for f in absrel_files:
        og = os.path.basename(f).split(".")[0]
        ogs.add(og)
        res, err = parse_absrel(f, test_branches)
        if res:
            absrel_by_og[og] = res
        else:
            absrel_errors.append((og, err))

    with open(out_path, "w") as o:
        o.write("OG\tRELAX_K\tRELAX_LR\tRELAX_p\tABSREL_min_p\tABSREL_min_p_branch\t"
                "ABSREL_n_test_branches_p05\tABSREL_n_test_branches_total\n")
        for og in sorted(ogs):
            r = relax_by_og.get(og, {})
            a = absrel_by_og.get(og, {})
            o.write(f"{og}\t"
                    f"{r.get('K', 'NA')}\t{r.get('LR', 'NA')}\t{r.get('p_value', 'NA')}\t"
                    f"{a.get('min_p', 'NA')}\t{a.get('min_p_branch', 'NA')}\t"
                    f"{a.get('n_test_branches_p05', 'NA')}\t{a.get('n_test_branches_total', 'NA')}\n")

    print(f"\nWrote {out_path}: {len(ogs)} OGs "
          f"({len(relax_by_og)} with parsed RELAX, {len(absrel_by_og)} with parsed aBSREL)")

    if relax_errors:
        print(f"\n{len(relax_errors)} RELAX JSONs failed to parse (first 5):")
        for og, err in relax_errors[:5]:
            print(f"  {og}: {err}")
        print("If this is most/all of them, the JSON schema doesn't match what this "
              "script expects -- send me one raw .RELAX.json and I'll fix it.")
    if absrel_errors:
        print(f"\n{len(absrel_errors)} aBSREL JSONs failed to parse (first 5):")
        for og, err in absrel_errors[:5]:
            print(f"  {og}: {err}")
        print("If this is most/all of them, send me one raw .ABSREL.json and I'll fix it.")


if __name__ == "__main__":
    main()
