#!/usr/bin/env bash
set -euo pipefail

# codeml_batch_launcher.sh
#
# Creates one directory per OG, runs codeml null (one-ratio) and alt (two-ratio branch model),
# parses lnL, computes LRT (df=1), and writes a summary TSV.
#
# Requirements:
# - codeml in PATH
# - A foreground-labeled Newick tree for PAML with HetGla and HelKap terminal branches marked #1
#
# Inputs:
#   1) PAML alignment directory containing per-OG PHYLIP codon alignments: paml_phy/OGxxxx.codon.phy
#   2) Foreground-labeled PAML tree: species.paml_fg.nwk
#   3) Output root directory (will create one folder per OG)
#
# Example:
#   bash codeml_batch_launcher.sh /fast/.../selection_from_SCO/paml_phy species.paml_fg.nwk /fast/.../paml_runs

PAML_PHY_DIR="${1:?Usage: bash codeml_batch_launcher.sh <paml_phy_dir> <paml_fg_tree.nwk> <out_root>}"
TREE_FG="${2:?Usage: bash codeml_batch_launcher.sh <paml_phy_dir> <paml_fg_tree.nwk> <out_root>}"
OUTROOT="${3:?Usage: bash codeml_batch_launcher.sh <paml_phy_dir> <paml_fg_tree.nwk> <out_root>}"

CODEML="${CODEML:-codeml}"

command -v "${CODEML}" >/dev/null 2>&1 || { echo "codeml not found in PATH (CODEML=${CODEML})"; exit 2; }
[ -s "${TREE_FG}" ] || { echo "Missing tree: ${TREE_FG}"; exit 2; }

mkdir -p "${OUTROOT}"

SUMMARY="${OUTROOT}/codeml_branch_models_summary.tsv"
printf "OG\tntaxa\tnsites\tlnL_null\tlnL_alt\tLRT\tp_df1\tnull_out\talt_out\n" > "${SUMMARY}"

parse_phy_header() {
  local phy="$1"
  awk 'NR==1{print $1"\t"$2; exit}' "$phy"
}

parse_lnl() {
  local out="$1"
  awk '
    /lnL\(/ && /=/ {
      for(i=1;i<=NF;i++){
        if($i=="=" && (i+1)<=NF){
          print $(i+1);
          exit
        }
      }
    }
  ' "$out" | head -n1
}

pval_df1() {
  local lrt="$1"
  python3 - <<PY
import math
x=float("${lrt}")
p = math.erfc(math.sqrt(x/2.0)) if x >= 0 else 1.0
print(p)
PY
}

shopt -s nullglob
PHY_FILES=("${PAML_PHY_DIR}"/*.phy)
if [ "${#PHY_FILES[@]}" -eq 0 ]; then
  echo "No .phy files found in ${PAML_PHY_DIR}"
  exit 1
fi

for phy in "${PHY_FILES[@]}"; do
  base="$(basename "$phy")"
  og="${base%%.*}"
  ogdir="${OUTROOT}/${og}"
  mkdir -p "${ogdir}"

  cp -f "${phy}" "${ogdir}/alignment.phy"
  cp -f "${TREE_FG}" "${ogdir}/tree_fg.nwk"

  read -r ntaxa nsites < <(parse_phy_header "${ogdir}/alignment.phy")

  cat > "${ogdir}/codeml_null.ctl" <<'CTL'
      seqfile = alignment.phy
     treefile = tree_fg.nwk
      outfile = codeml_null.out

        noisy = 3
      verbose = 0
      runmode = 0

      seqtype = 1       * 1:codons
    CodonFreq = 2       * 2:F3x4

        clock = 0
       aaDist = 0

      model = 0         * 0: one-ratio
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

      seqtype = 1       * 1:codons
    CodonFreq = 2       * 2:F3x4

        clock = 0
       aaDist = 0

      model = 2         * 2: two-ratio branch model (foreground labeled #1)
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
    p="$(pval_df1 "${lrt}")"
  else
    lrt="NA"
    p="NA"
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"     "${og}" "${ntaxa}" "${nsites}"     "${lnl_null:-NA}" "${lnl_alt:-NA}" "${lrt}" "${p}"     "${ogdir}/codeml_null.out" "${ogdir}/codeml_alt.out" >> "${SUMMARY}"
done

echo "Wrote summary: ${SUMMARY}"
