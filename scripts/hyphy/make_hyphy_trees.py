#!/usr/bin/env python3
"""
Converts an existing PAML-style foreground tree (#1-tagged, as used by
codeml) into HyPhy-labeled trees for RELAX and aBSREL, using the SAME
foreground branches as the PAML run -- so RELAX/aBSREL test the identical
hypothesis as the branch-model codeml test on that phenotype, and results
are directly comparable.

Produces two output trees per input:
  <prefix>.RELAX.nwk      : every leaf labeled {Test} (the #1-tagged ones)
                             or {Reference} (everything else) -- every leaf
                             is labeled explicitly (never left unlabeled),
                             so the result doesn't depend on how a given
                             HyPhy version treats unlabeled branches.
  <prefix>.aBSREL_fg.nwk  : only the #1-tagged (foreground) leaves labeled
                             {Test}, everything else left unlabeled -- for
                             an aBSREL run restricted to just the phenotype
                             branches (--branches Test), which tests
                             specifically whether those branches show
                             episodic diversifying selection, instead of
                             paying the multiple-testing cost of testing
                             all 16 branches.

(A genome-wide, all-branches aBSREL run needs no relabeling at all --
just strip the #1 tags and use the plain topology; aBSREL tests every
branch by default. Use --branches ALL, or the default with no --branches
flag, against the plain tree if you also want that.)

Usage:
    python3 make_hyphy_trees.py <phenotype_paml_fg.nwk> <output_prefix>

e.g.
    python3 make_hyphy_trees.py acid_species.paml_fg.nwk hyphy_trees/acid
      -> hyphy_trees/acid.RELAX.nwk
      -> hyphy_trees/acid.aBSREL_fg.nwk
"""
import sys
import os
import re


def parse_leaves_and_tags(newick_text):
    """Return the raw text with every '<name>#1' or '<name>' leaf token
    identified, without touching internal-node labels/support values."""
    # leaf tokens: alnum/underscore name, optionally #1, immediately
    # followed by ':' (branch length) -- same convention as validate_trees.py
    return re.findall(r"([A-Za-z][A-Za-z0-9_]*)(#1)?(?=:)", newick_text)


def build_relax_tree(text):
    # NB: the (?=:) is a lookahead -- it does not consume the ':', so the
    # replacement must NOT append one itself (that was bug: produced "::").
    def repl(m):
        name, tag = m.group(1), m.group(2)
        label = "Test" if tag else "Reference"
        return f"{name}{{{label}}}"
    return re.sub(r"([A-Za-z][A-Za-z0-9_]*)(#1)?(?=:)", repl, text)


def build_absrel_fg_tree(text):
    def repl(m):
        name, tag = m.group(1), m.group(2)
        if tag:
            return f"{name}{{Test}}"
        return name
    return re.sub(r"([A-Za-z][A-Za-z0-9_]*)(#1)?(?=:)", repl, text)


def main():
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <phenotype_paml_fg.nwk> <output_prefix>")
    in_path, out_prefix = sys.argv[1:3]

    text = open(in_path).read().strip()
    # strip a leading tree-count header line if present (phylip-style)
    lines = [l for l in text.splitlines() if l.strip() and not re.fullmatch(r"\d+", l.strip())]
    text = "\n".join(lines)

    leaves = parse_leaves_and_tags(text)
    fg = [name for name, tag in leaves if tag]
    bg = [name for name, tag in leaves if not tag]
    if not fg:
        sys.exit(f"No #1-tagged foreground leaves found in {in_path} -- nothing to convert.")
    print(f"Input: {in_path}")
    print(f"  Foreground (#1) leaves -> Test:      {fg}")
    print(f"  Background leaves      -> Reference: {bg}")

    relax_tree = build_relax_tree(text)
    absrel_tree = build_absrel_fg_tree(text)

    relax_path = f"{out_prefix}.RELAX.nwk"
    absrel_path = f"{out_prefix}.aBSREL_fg.nwk"
    out_dir = os.path.dirname(relax_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(relax_path, "w") as o:
        o.write(relax_tree + "\n")
    with open(absrel_path, "w") as o:
        o.write(absrel_tree + "\n")

    print(f"\nWrote {relax_path}  (RELAX: --test Test --reference Reference)")
    print(f"Wrote {absrel_path}  (aBSREL: --branches Test)")


if __name__ == "__main__":
    main()
