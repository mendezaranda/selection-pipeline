#!/usr/bin/env python3
"""
Prune a HyPhy-labeled tree (leaves tagged {Test}/{Reference} by
make_hyphy_trees.py) down to exactly the taxa present in one OG's codon
alignment, preserving each leaf's {Label} tag.

Same problem, same fix as bug 3 in og_to_codon_paml_codeml.sh: not every OG
has all 16 species (missing single-copy orthologs is normal), and HyPhy
requires an exact match between tree tips and alignment sequences -- reusing
one static full-taxa tree for every OG makes HyPhy abort with "The number
of tree tips ... is not equal to the number of sequences ..." for any OG
missing even one species. This prunes per-OG instead, the same way codeml's
tree gets pruned per-OG via TREE_FG.

Usage:
    python3 prune_hyphy_tree.py <codon_fasta> <tree_in.nwk> <tree_out.nwk> [required_tags]

<codon_fasta>   : OG's hyphy_fa/<OG>.codon.fasta -- its headers (species
                  codes) are the taxon set to prune down to.
<tree_in.nwk>   : full-taxa RELAX or aBSREL_fg tree from make_hyphy_trees.py
<tree_out.nwk>  : pruned tree, written only on success
[required_tags] : comma-separated tag names that must have >=1 leaf left
                  after pruning, else treated as a skip, not a hard error.
                  Default "Test". Pass "Test,Reference" for a RELAX tree
                  (RELAX needs both categories present).

Exit codes: 0 OK, 2 a required tag has zero leaves left after pruning
(e.g. this OG's alignment is missing every Test-branch species), 3 no taxa
left at all, 4 internal taxon-set mismatch after pruning (should not
happen; would indicate a bug in this script).
"""
import sys
import re


def read_fasta_names(fn):
    names = []
    for line in open(fn):
        if line.startswith(">"):
            names.append(line[1:].strip().split()[0])
    return names


def parse_newick(s):
    s = s.strip()
    if s.endswith(';'):
        s = s[:-1]
    pos = [0]

    def parse_label():
        start = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ',();:':
            pos[0] += 1
        return s[start:pos[0]]

    def split_tag(label):
        # leaf/internal labels here look like "HelKap{Test}" or "HelKap"
        # (never "#1" -- that's PAML's convention, not HyPhy's)
        m = re.match(r'^(.*?)(\{[A-Za-z_][A-Za-z0-9_]*\})?$', label)
        return m.group(1), m.group(2)

    def parse_node():
        node = {'name': None, 'tag': None, 'length': None, 'children': []}
        if s[pos[0]] == '(':
            pos[0] += 1
            node['children'].append(parse_node())
            while s[pos[0]] == ',':
                pos[0] += 1
                node['children'].append(parse_node())
            assert s[pos[0]] == ')', f"expected ) at {pos[0]}"
            pos[0] += 1
            label = parse_label()
            name, tag = split_tag(label)
            node['name'] = name if name else None
            node['tag'] = tag
        else:
            label = parse_label()
            name, tag = split_tag(label)
            node['name'] = name
            node['tag'] = tag
        if pos[0] < len(s) and s[pos[0]] == ':':
            pos[0] += 1
            node['length'] = parse_label()
        return node

    return parse_node()


def leaves(node, out=None):
    if out is None:
        out = []
    if not node['children']:
        out.append(node)
    else:
        for c in node['children']:
            leaves(c, out)
    return out


def prune(node, keep):
    if not node['children']:
        return dict(node) if node['name'] in keep else None
    new_children = []
    for c in node['children']:
        pc = prune(c, keep)
        if pc is not None:
            new_children.append(pc)
    if not new_children:
        return None
    if len(new_children) == 1:
        # Collapse a single-child internal node into its child, merging
        # branch lengths. The child's own {Label} tag (if it's a leaf) is
        # already on the child dict, so it's preserved automatically --
        # unlike codeml's #1 pruning, there's no separate "transfer the
        # foreground tag" step because every leaf already carries its own
        # tag going in.
        child = new_children[0]

        def add_len(a, b):
            fa = float(a) if a not in (None, '') else 0.0
            fb = float(b) if b not in (None, '') else 0.0
            return repr(fa + fb)
        merged_len = add_len(node.get('length'), child.get('length'))
        child = dict(child)
        child['length'] = merged_len
        return child
    return {'name': node['name'], 'tag': node['tag'], 'length': node['length'], 'children': new_children}


def to_newick(node):
    def tag_str(n):
        return n['tag'] if n['tag'] else ''

    def rec(n):
        if not n['children']:
            s = (n['name'] or '') + tag_str(n)
        else:
            s = '(' + ','.join(rec(c) for c in n['children']) + ')' + (n['name'] or '') + tag_str(n)
        if n.get('length') not in (None, ''):
            s += ':' + n['length']
        return s
    return rec(node) + ';'


def main():
    if len(sys.argv) not in (4, 5):
        sys.exit(f"usage: {sys.argv[0]} <codon_fasta> <tree_in.nwk> <tree_out.nwk> [required_tags]")
    codon_fa, tree_in, tree_out = sys.argv[1:4]
    required_tags = sys.argv[4].split(",") if len(sys.argv) == 5 else ["Test"]

    keep = set(read_fasta_names(codon_fa))

    tree = parse_newick(open(tree_in).read())
    pruned = prune(tree, keep)
    if pruned is None:
        print("NO_TAXA_LEFT")
        sys.exit(3)

    pruned_leaves = leaves(pruned)
    pruned_names = set(l['name'] for l in pruned_leaves)
    if pruned_names != keep:
        print(f"TAXON_SET_MISMATCH: pruned={sorted(pruned_names)} keep={sorted(keep)}")
        sys.exit(4)

    for tag in required_tags:
        tagged = [l for l in pruned_leaves if l['tag'] == "{" + tag + "}"]
        if not tagged:
            print(f"NO_{tag.upper()}_TAXA_LEFT")
            sys.exit(2)

    with open(tree_out, "w") as o:
        o.write(to_newick(pruned) + "\n")
    counts = ", ".join(
        f"{tag}={sum(1 for l in pruned_leaves if l['tag'] == '{' + tag + '}')}"
        for tag in required_tags
    )
    print(f"OK: pruned to {len(pruned_names)} taxa ({counts})")


if __name__ == "__main__":
    main()
