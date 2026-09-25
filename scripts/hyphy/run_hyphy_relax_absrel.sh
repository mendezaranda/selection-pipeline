#!/usr/bin/env bash
set -eo pipefail
#
# Runs HyPhy RELAX and/or aBSREL (branches-restricted to the phenotype
# foreground) for every OG a phenotype's codeml run already produced, reusing
# the existing hyphy_fa/<OG>.codon.fasta files -- these are byte-identical to
# codon_aln (og_to_codon_paml_codeml.sh just copies codon_fa -> hyphy_fa), so
# no realignment, no rerun of the expensive alignment step.
#
# Both tests use the SAME Test/foreground branches as the phenotype's
# existing codeml branch-model tree (built by make_hyphy_trees.py from that
# tree), so results are directly comparable to the codeml LRT for the same
# phenotype:
#   RELAX  : --test Test --reference Reference  (all 16 taxa split into the
#            two sets -- tests whether selection strength is relaxed [k<1]
#            or intensified [k>1] on the phenotype branches vs background)
#   aBSREL : --branches Test  (tests only the phenotype branches for episodic
#            diversifying selection, instead of all 16 branches -- keeps the
#            multiple-testing burden the same as testing the phenotype
#            hypothesis specifically, not a genome-wide all-branch scan; also,
#            unlike codeml's two-ratio model, aBSREL fits each Test branch's
#            own rate distribution independently -- it does NOT pool omega
#            across the Test set, so it doesn't have the dilution issue that
#            makes codeml's pooled-foreground results hard to attribute to a
#            single lineage)
#
# PER-OG TREE PRUNING (added 2026-09-25, bug fix): not every OG has all 16
# species (missing single-copy orthologs is normal -- this is the exact same
# problem bug 3 in og_to_codon_paml_codeml.sh fixed for codeml). HyPhy
# requires an EXACT match between tree tips and alignment sequences, so
# reusing RELAX_TREE/ABSREL_TREE as-is for every OG made HyPhy abort ("The
# number of tree tips ... is not equal to the number of sequences ...") for
# any OG missing even one species relative to the full-taxa tree -- silently,
# per-OG, leaving an empty output JSON and only a message in that OG's log.
# Fixed by pruning each tree down to the OG's actual taxon set (via
# prune_hyphy_tree.py, same pure-Python Newick pruner approach as codeml's
# fix) before every hyphy call. An OG whose alignment has lost every Test (or,
# for RELAX, every Reference) taxon is logged and skipped for that test
# rather than passed to HyPhy, which would just error on it anyway.
#
# ONLY=absrel|relax|both (env var, default: both) skips the compute (not
# just the reporting) for whichever method you don't want yet -- e.g. run
# aBSREL first, add RELAX later with the same command + ONLY=relax,
# already-done aBSREL JSONs are left alone (skip-if-exists, like the rest of
# this pipeline). This is an env var, not a --flag, on purpose: a leading
# "--" in a copy-pasted command can get silently mangled into an en-dash by
# some renderers, which breaks flag matching with no obvious error -- an env
# var assignment can't be corrupted that way.
#
# Usage:
#   ONLY=absrel bash run_hyphy_relax_absrel.sh \
#       <PHENOTYPE_DIR> <RELAX_TREE.nwk> <ABSREL_FG_TREE.nwk> <OUTDIR> [OG_LIST]
#
# <PHENOTYPE_DIR>     : e.g. 06_paml/acid  (needs hyphy_fa/<OG>.codon.fasta)
# <RELAX_TREE.nwk>    : from make_hyphy_trees.py, e.g. acid.RELAX.nwk
#                        (still required even with ONLY=absrel -- just used
#                        to skip the RELAX call, not to avoid generating it)
# <ABSREL_FG_TREE.nwk>: from make_hyphy_trees.py, e.g. acid.aBSREL_fg.nwk
# <OUTDIR>            : where results go, e.g. 06_paml/acid/hyphy_runs
# [OG_LIST]           : optional -- one OG id per line, to restrict to a
#                        subset (e.g. only the OGs significant in codeml).
#                        Default: every hyphy_fa/*.codon.fasta in
#                        PHENOTYPE_DIR.
#
# HYPHY executable: set via the HYPHY env var if it isn't just "hyphy" on
# your PATH, e.g.:
#   export HYPHY=/home/YOUR_USERNAME/hyphy/bin/hyphy
# (same override pattern as CODEML in codeml_from_existing_alignment.sh --
# I don't know your cluster's HyPhy install path/module, so check
# `which hyphy` or `module avail hyphy` first and set this if needed.)

ONLY="${ONLY:-both}"
case "${ONLY}" in
  absrel|relax|both) ;;
  *) echo "ONLY must be absrel, relax, or both (got '${ONLY}')"; exit 2 ;;
esac

PHENOTYPE_DIR="${1:?Usage: [ONLY=absrel|relax|both] bash run_hyphy_relax_absrel.sh <PHENOTYPE_DIR> <RELAX_TREE.nwk> <ABSREL_FG_TREE.nwk> <OUTDIR> [OG_LIST]}"
RELAX_TREE="${2:?Usage: [ONLY=absrel|relax|both] bash run_hyphy_relax_absrel.sh <PHENOTYPE_DIR> <RELAX_TREE.nwk> <ABSREL_FG_TREE.nwk> <OUTDIR> [OG_LIST]}"
ABSREL_TREE="${3:?Usage: [ONLY=absrel|relax|both] bash run_hyphy_relax_absrel.sh <PHENOTYPE_DIR> <RELAX_TREE.nwk> <ABSREL_FG_TREE.nwk> <OUTDIR> [OG_LIST]}"
OUTDIR="${4:?Usage: [ONLY=absrel|relax|both] bash run_hyphy_relax_absrel.sh <PHENOTYPE_DIR> <RELAX_TREE.nwk> <ABSREL_FG_TREE.nwk> <OUTDIR> [OG_LIST]}"
OG_LIST="${5:-}"

HYPHY="${HYPHY:-hyphy}"
command -v "${HYPHY}" >/dev/null 2>&1 || {
  echo "'${HYPHY}' not found on PATH. Set HYPHY=/path/to/hyphy (or 'module load hyphy' first) and rerun."
  exit 2
}

PRUNE_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/prune_hyphy_tree.py"
[ -s "${PRUNE_SCRIPT}" ] || { echo "Missing prune_hyphy_tree.py (expected next to this script at ${PRUNE_SCRIPT})"; exit 2; }

HYPHY_FA_DIR="${PHENOTYPE_DIR}/hyphy_fa"
[ -d "${HYPHY_FA_DIR}" ] || { echo "Missing ${HYPHY_FA_DIR}"; exit 2; }
[ -s "${RELAX_TREE}" ] || { echo "Missing RELAX tree: ${RELAX_TREE}"; exit 2; }
[ -s "${ABSREL_TREE}" ] || { echo "Missing aBSREL tree: ${ABSREL_TREE}"; exit 2; }

mkdir -p "${OUTDIR}/relax" "${OUTDIR}/absrel" "${OUTDIR}/logs" "${OUTDIR}/trees"

if [ -n "${OG_LIST}" ]; then
  mapfile -t OGS < "${OG_LIST}"
else
  OGS=()
  for f in "${HYPHY_FA_DIR}"/*.codon.fasta; do
    [ -e "${f}" ] || continue
    base="$(basename "${f}")"
    OGS+=("${base%%.*}")
  done
fi
echo "Running HyPhy (ONLY=${ONLY}) for ${#OGS[@]} OGs from ${HYPHY_FA_DIR}"

n_ok=0
n_skip=0
n_pruneskip=0
for og in "${OGS[@]}"; do
  [ -z "${og}" ] && continue
  fa="${HYPHY_FA_DIR}/${og}.codon.fasta"
  if [ ! -s "${fa}" ]; then
    echo "SKIP ${og}: no ${fa}"
    n_skip=$((n_skip + 1))
    continue
  fi

  if [ "${ONLY}" = "relax" ] || [ "${ONLY}" = "both" ]; then
    relax_out="${OUTDIR}/relax/${og}.RELAX.json"
    if [ ! -s "${relax_out}" ]; then
      relax_tree_pruned="${OUTDIR}/trees/${og}.RELAX.pruned.nwk"
      prune_log="${OUTDIR}/logs/${og}.relax.prune.log"
      if python3 "${PRUNE_SCRIPT}" "${fa}" "${RELAX_TREE}" "${relax_tree_pruned}" "Test,Reference" > "${prune_log}" 2>&1; then
        "${HYPHY}" relax --alignment "${fa}" --tree "${relax_tree_pruned}" \
          --test Test --reference Reference \
          --output "${relax_out}" \
          > "${OUTDIR}/logs/${og}.relax.log" 2>&1 || echo "  RELAX FAILED: ${og} (see ${OUTDIR}/logs/${og}.relax.log)"
      else
        echo "  RELAX SKIP: ${og} (tree pruning: $(cat "${prune_log}"))"
        n_pruneskip=$((n_pruneskip + 1))
      fi
    fi
  fi

  if [ "${ONLY}" = "absrel" ] || [ "${ONLY}" = "both" ]; then
    absrel_out="${OUTDIR}/absrel/${og}.ABSREL.json"
    if [ ! -s "${absrel_out}" ]; then
      absrel_tree_pruned="${OUTDIR}/trees/${og}.aBSREL_fg.pruned.nwk"
      prune_log="${OUTDIR}/logs/${og}.absrel.prune.log"
      if python3 "${PRUNE_SCRIPT}" "${fa}" "${ABSREL_TREE}" "${absrel_tree_pruned}" "Test" > "${prune_log}" 2>&1; then
        "${HYPHY}" absrel --alignment "${fa}" --tree "${absrel_tree_pruned}" \
          --branches Test \
          --output "${absrel_out}" \
          > "${OUTDIR}/logs/${og}.absrel.log" 2>&1 || echo "  aBSREL FAILED: ${og} (see ${OUTDIR}/logs/${og}.absrel.log)"
      else
        echo "  aBSREL SKIP: ${og} (tree pruning: $(cat "${prune_log}"))"
        n_pruneskip=$((n_pruneskip + 1))
      fi
    fi
  fi

  n_ok=$((n_ok + 1))
done

echo
echo "=== Done: ${n_ok} OGs attempted, ${n_skip} skipped (no fasta), ${n_pruneskip} test(s) skipped (tree pruning: OG missing every taxon on one side) ==="
[ "${ONLY}" = "relax" ] || [ "${ONLY}" = "both" ] && echo "RELAX JSONs:  ${OUTDIR}/relax/<OG>.RELAX.json"
[ "${ONLY}" = "absrel" ] || [ "${ONLY}" = "both" ] && echo "aBSREL JSONs: ${OUTDIR}/absrel/<OG>.ABSREL.json"
echo "Per-OG pruned trees: ${OUTDIR}/trees/<OG>.{RELAX,aBSREL_fg}.pruned.nwk"
echo "Next: parse_hyphy_results.py will pull p-values / K (relaxation param) out of these into one TSV."
