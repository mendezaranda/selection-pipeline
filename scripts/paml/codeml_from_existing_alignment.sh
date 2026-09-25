#!/usr/bin/env bash
set -eo pipefail

# Rerun ONLY the tree-pruning + codeml branch-model step for one OG, reusing an
# already-built codon alignment (paml_phy) from a prior full run (e.g.
# only_HelKap, built by build_codon_alignments_array.sh). This skips
# mafft/trimAl/backtranslation entirely, since those don't depend on which
# taxa are marked foreground -- only the tree does.
#
# Usage:
#   bash codeml_from_existing_alignment.sh <PAML_PHY> <TREE_FG.nwk> <OUTDIR> [THREADS]
#
# <PAML_PHY> is an existing .codon.phy file, e.g.
#   06_paml/only_HelKap/paml_phy/OG0007353.codon.phy
# from a previous og_to_codon_paml_codeml.sh run. Its taxon set is whatever
# that run's MIN_TAXA/filtering already decided for this OG -- this script
# does not re-filter or re-align anything.

PAML_PHY_SRC="${1:?Usage: bash codeml_from_existing_alignment.sh <PAML_PHY> <TREE_FG.nwk> <OUTDIR> [THREADS]}"
TREE_FG="${2:?Usage: bash codeml_from_existing_alignment.sh <PAML_PHY> <TREE_FG.nwk> <OUTDIR> [THREADS]}"
OUTDIR="${3:?Usage: bash codeml_from_existing_alignment.sh <PAML_PHY> <TREE_FG.nwk> <OUTDIR> [THREADS]}"

CODEML="${CODEML:-codeml}"
command -v "${CODEML}" >/dev/null 2>&1 || { echo "codeml not found in PATH (CODEML=${CODEML})"; exit 2; }
[ -s "${TREE_FG}" ] || { echo "Missing tree file: ${TREE_FG}"; exit 2; }
[ -s "${PAML_PHY_SRC}" ] || { echo "Missing source alignment: ${PAML_PHY_SRC}"; exit 2; }

mkdir -p "${OUTDIR}/codeml_runs" "${OUTDIR}/paml_phy_used"

ogbase="$(basename "${PAML_PHY_SRC}")"
ogname="${ogbase%%.*}"
ogdir="${OUTDIR}/codeml_runs/${ogname}"

TASK_TAG="${SLURM_ARRAY_TASK_ID:-0}"
SUMMARY_TSV="${OUTDIR}/codeml_summary.task${TASK_TAG}.tsv"

append_summary () {
  local og="$1" st="$2" ntaxa="$3" nsites="$4" ln0="$5" ln1="$6" lrt="$7" p="$8" note="$9"
  if [ ! -s "${SUMMARY_TSV}" ]; then
    echo -e "OG\tstatus\tntaxa\tnsites\tlnL_null\tlnL_alt\tLRT\tp_df1\tnote" > "${SUMMARY_TSV}"
  fi
  echo -e "${og}\t${st}\t${ntaxa}\t${nsites}\t${ln0}\t${ln1}\t${lrt}\t${p}\t${note}" >> "${SUMMARY_TSV}"
}

# Same lnL parser as og_to_codon_paml_codeml.sh (see that script for why this
# format, not a literal "=", is required).
parse_lnl() {
  local f="$1" v
  v="$(grep -m1 -oE 'lnL\([^)]*\):[[:space:]]*-?[0-9]+\.[0-9]+' "$f" 2>/dev/null \
       | grep -oE -- '-?[0-9]+\.[0-9]+$')"
  if [ -z "${v}" ]; then
    v="$(awk '/lnL/ && /=/{for(i=1;i<=NF;i++){if($i=="="){print $(i+1); exit}}}' "$f" 2>/dev/null | head -n1)"
  fi
  echo "${v}"
}

codeml_alt_out="${ogdir}/codeml_alt.out"
codeml_ok=0
if [ -s "${codeml_alt_out}" ] && [ -n "$(parse_lnl "${codeml_alt_out}")" ]; then
  codeml_ok=1
fi
if [ "${codeml_ok}" -eq 1 ] || [ -e "${ogdir}/SKIPPED_NO_FOREGROUND" ]; then
  echo "SKIP ${ogname} (exists)"
  exit 0
fi

mkdir -p "${ogdir}"
cp -f "${PAML_PHY_SRC}" "${ogdir}/alignment.phy"
read -r ntaxa nsites < <(awk 'NR==1{print $1"\t"$2; exit}' "${ogdir}/alignment.phy")

# Species present in this OG = the taxon names already baked into the phylip
# alignment (first column of each sequence line, after the 1-line header).
species_list="${ogdir}/species.txt"
awk 'NR>1{print $1}' "${ogdir}/alignment.phy" > "${species_list}"

# Prune TREE_FG down to exactly these taxa, same logic as bug 3's fix in
# og_to_codon_paml_codeml.sh (kept identical so both scripts stay consistent).
if python3 - "${species_list}" "${TREE_FG}" "${ogdir}/tree_fg.nwk" <<'PY' > "${ogdir}/prune_tree.log" 2>&1
import sys, re

SPECIES_FILE, TREE_IN, TREE_OUT = sys.argv[1], sys.argv[2], sys.argv[3]
keep = set(x.strip() for x in open(SPECIES_FILE) if x.strip())

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
        m = re.match(r'^(.*?)(#\d+)?$', label)
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
    if out is None: out = []
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

tree = parse_newick(open(TREE_IN).read())
pruned = prune(tree, keep)
if pruned is None:
    print("NO_TAXA_LEFT")
    sys.exit(3)

pruned_leaves = leaves(pruned)
pruned_names = set(l['name'] for l in pruned_leaves)
if pruned_names != keep:
    print(f"TAXON_SET_MISMATCH: pruned={sorted(pruned_names)} keep={sorted(keep)}")
    sys.exit(4)

fg = [l for l in pruned_leaves if l['tag']]
if not fg:
    print("NO_FOREGROUND_TAXON")
    sys.exit(2)

with open(TREE_OUT, "w") as o:
    o.write(to_newick(pruned) + "\n")
print(f"OK: pruned to {len(pruned_names)} taxa, foreground={[l['name'] for l in fg]}")
PY
then
  prune_rc=0
else
  prune_rc=$?
fi

if [ "${prune_rc}" -eq 2 ]; then
  echo "SKIP ${ogname} (foreground taxon absent from this OG's alignment; see ${ogdir}/prune_tree.log)"
  touch "${ogdir}/SKIPPED_NO_FOREGROUND"
  append_summary "${ogname}" "SKIP_NO_FOREGROUND_TAXON" "${ntaxa}" "${nsites}" "NA" "NA" "NA" "NA" "foreground_taxon_missing"
  exit 0
elif [ "${prune_rc}" -ne 0 ]; then
  echo "SKIP ${ogname} (tree pruning failed rc=${prune_rc}; see ${ogdir}/prune_tree.log)"
  append_summary "${ogname}" "SKIP_TREE_PRUNE_FAILED" "${ntaxa}" "${nsites}" "NA" "NA" "NA" "NA" "tree_prune_failed"
  exit 0
fi

cat > "${ogdir}/codeml_null.ctl" <<'CTL'
      seqfile = alignment.phy
     treefile = tree_fg.nwk
      outfile = codeml_null.out

        noisy = 3
      verbose = 0
      runmode = 0

      seqtype = 1
    CodonFreq = 2

        clock = 0
       aaDist = 0

      model = 0
    NSsites = 0

  fix_kappa = 0
      kappa = 2

  fix_omega = 0
      omega = 0.2

    cleandata = 0
CTL

cat > "${ogdir}/codeml_alt.ctl" <<'CTL'
      seqfile = alignment.phy
     treefile = tree_fg.nwk
      outfile = codeml_alt.out

        noisy = 3
      verbose = 0
      runmode = 0

      seqtype = 1
    CodonFreq = 2

        clock = 0
       aaDist = 0

      model = 2
    NSsites = 0

  fix_kappa = 0
      kappa = 2

  fix_omega = 0
      omega = 0.2

    cleandata = 0
CTL

( cd "${ogdir}" && "${CODEML}" codeml_null.ctl > codeml_null.log 2>&1 ) || true
( cd "${ogdir}" && "${CODEML}" codeml_alt.ctl  > codeml_alt.log  2>&1 ) || true

lnl_null="$(parse_lnl "${ogdir}/codeml_null.out" || true)"
lnl_alt="$(parse_lnl "${ogdir}/codeml_alt.out"  || true)"

if [ -n "${lnl_null}" ] && [ -n "${lnl_alt}" ]; then
  lrt="$(python3 - <<PY
ln0=float("${lnl_null}")
ln1=float("${lnl_alt}")
print(max(0.0, 2.0*(ln1-ln0)))
PY
)"
  p="$(python3 - <<PY
import math
x=float("${lrt}")
p = math.erfc(math.sqrt(x/2.0)) if x >= 0 else 1.0
print(p)
PY
)"
  append_summary "${ogname}" "OK" "${ntaxa}" "${nsites}" "${lnl_null}" "${lnl_alt}" "${lrt}" "${p}" "codeml_branch_model"
else
  append_summary "${ogname}" "CODEML_FAIL" "${ntaxa}" "${nsites}" "${lnl_null:-NA}" "${lnl_alt:-NA}" "NA" "NA" "parse_lnl_failed"
fi

echo "OK ${ogname} -> ${ogdir}"
