#!/usr/bin/env bash
#SBATCH --job-name=codeml_reuse
#SBATCH -o /home/%u/outdir/codeml_reuse_%A_%a.out
#SBATCH -e /home/%u/errdir/codeml_reuse_%A_%a.err
#SBATCH -t 08:00:00
#SBATCH -N1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --array=1-20

set -eo pipefail

# Reruns ONLY the tree-pruning + codeml branch-model step for a new
# phenotype's foreground tree, reusing codon alignments already built by a
# prior full build_codon_alignments_array.sh run (default source:
# only_HelKap). No mafft/trimAl, so this needs far less time/memory per task
# than the full pipeline -- adjust the #SBATCH lines above if your
# alignments are unusually large.
#
# source configs/env.sh (see configs/env.example.sh) first, or export
# RUN_BASE yourself. Then set PHENOTYPE and TREE_FG (or override via
# --export), e.g.:
#   sbatch --export=ALL,PHENOTYPE=acid,TREE_FG=/path/to/acid/species.paml_fg.nwk \
#       run_codeml_from_existing_array.sh

PHENOTYPE="${PHENOTYPE:?Set PHENOTYPE, e.g. acid | acid_wBat | capsaicin | HelKap_HetGla}"
RUN_BASE="${RUN_BASE:?Set RUN_BASE (see configs/env.example.sh)}"
PAML_DIR="${PAML_DIR:-${RUN_BASE}/06_paml}"

SRC_PAML_PHY_DIR="${SRC_PAML_PHY_DIR:-${PAML_DIR}/only_HelKap/paml_phy}"
TREE_FG="${TREE_FG:?Set TREE_FG to the phenotype-specific foreground tree, e.g. ${PAML_DIR}/${PHENOTYPE}/species.paml_fg.nwk}"
OUTDIR="${OUTDIR:-${PAML_DIR}/${PHENOTYPE}}"

# Self-locating, like every other array driver in this repo -- works
# wherever scripts/paml/ is cloned or copied to.
SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER="${WORKER:-${SCRIPTS_DIR}/codeml_from_existing_alignment.sh}"
CODEML="${CODEML:-codeml}"

mkdir -p "${OUTDIR}"
[ -d "${SRC_PAML_PHY_DIR}" ] || { echo "Missing SRC_PAML_PHY_DIR: ${SRC_PAML_PHY_DIR}"; exit 1; }
[ -s "${TREE_FG}" ] || { echo "Missing TREE_FG: ${TREE_FG}"; exit 1; }
[ -x "${WORKER}" ] || { echo "Missing/not executable WORKER: ${WORKER}"; exit 1; }

mapfile -t PHYFILES < <(ls "${SRC_PAML_PHY_DIR}"/OG*.codon.phy 2>/dev/null | sort)
N=${#PHYFILES[@]}
[ "${N}" -gt 0 ] || { echo "No *.codon.phy found in ${SRC_PAML_PHY_DIR}"; exit 1; }

TASK_ID="${SLURM_ARRAY_TASK_ID}"
STEP="${SLURM_ARRAY_TASK_COUNT:-20}"
START=$((TASK_ID-1))

echo "PHENOTYPE=${PHENOTYPE}  TREE_FG=${TREE_FG}  OUTDIR=${OUTDIR}"
echo "Total alignments: ${N}"
echo "Task ${TASK_ID}: processing indices ${START}, ${START}+${STEP}, ..."

for ((i=START; i<N; i+=STEP)); do
  phy="${PHYFILES[$i]}"
  CODEML="${CODEML}" bash "${WORKER}" "${phy}" "${TREE_FG}" "${OUTDIR}"
done

# Merge per-task summaries (task 1 only -- see the caveat about this not
# actually waiting for other tasks in build_codon_alignments_array.sh; for
# a guaranteed-complete merge, run this as a separate --dependency=afterany
# job instead).
if [ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]; then
  merged="${OUTDIR}/codeml_branch_models_summary.tsv"
  first=1
  rm -f "$merged"
  for f in "${OUTDIR}"/codeml_summary.task*.tsv; do
    [ -s "$f" ] || continue
    if [ $first -eq 1 ]; then
      cat "$f" > "$merged"
      first=0
    else
      tail -n +2 "$f" >> "$merged"
    fi
  done
  echo "Merged summary: $merged"
fi
