#!/usr/bin/env bash
#SBATCH --job-name=hyphy_relax_absrel
#SBATCH -o /home/%u/outdir/hyphy_%A_%a.out
#SBATCH -e /home/%u/errdir/hyphy_%A_%a.err
#SBATCH -t 24:00:00
#SBATCH -N1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-40

set -eo pipefail

# SLURM array driver for run_hyphy_relax_absrel.sh: splits the full OG list
# for a phenotype into SLURM_ARRAY_TASK_COUNT interleaved chunks (index
# stride, not contiguous blocks -- keeps chunks even in size and load even
# if some OGs are much slower than others) and hands each chunk to
# run_hyphy_relax_absrel.sh as an OG_LIST subset. That script already skips
# any OG whose output JSON exists, so a rerun (more array tasks, a retry
# after a timeout, ONLY=relax after an earlier ONLY=absrel pass) is
# safe/cheap.
#
# Genome-wide, not just a significant subset -- run with no extra OG_LIST
# restriction, same OG universe as this phenotype's codeml run.
#
# source configs/env.sh (see configs/env.example.sh) first, or export
# RUN_BASE yourself. Then set PHENOTYPE (required) before submitting, e.g.:
#   sbatch --export=ALL,PHENOTYPE=only_HelKap run_hyphy_relax_absrel_array.sh
#   sbatch --export=ALL,PHENOTYPE=acid,ONLY=absrel run_hyphy_relax_absrel_array.sh
#
# Run once per phenotype (only_HelKap, HelKap_HetGla, capsaicin, acid,
# acid_wBat, ...) -- submit one array job per phenotype, same as how the
# PAML phenotype reruns are each their own sbatch.
#
# Prerequisite: make_hyphy_trees.py must already have been run for this
# phenotype (produces the two tree files this script checks for below).

PHENOTYPE="${PHENOTYPE:?Set PHENOTYPE, e.g. only_HelKap | HelKap_HetGla | capsaicin | acid | acid_wBat}"
RUN_BASE="${RUN_BASE:?Set RUN_BASE (see configs/env.example.sh)}"
PAML_DIR="${PAML_DIR:-${RUN_BASE}/06_paml}"
HYPHY_DIR="${HYPHY_DIR:-${RUN_BASE}/07_hyphy}"

PHENOTYPE_DIR="${PHENOTYPE_DIR:-${PAML_DIR}/${PHENOTYPE}}"
HYPHY_BASE="${HYPHY_BASE:-${HYPHY_DIR}/${PHENOTYPE}}"
RELAX_TREE="${RELAX_TREE:-${HYPHY_BASE}/hyphy_trees/${PHENOTYPE}.RELAX.nwk}"
ABSREL_TREE="${ABSREL_TREE:-${HYPHY_BASE}/hyphy_trees/${PHENOTYPE}.aBSREL_fg.nwk}"
OUTDIR="${OUTDIR:-${HYPHY_BASE}/hyphy_runs}"
TASKLIST_DIR="${HYPHY_BASE}/hyphy_runs/tasklists"

# Self-locating, like every other array driver in this repo.
SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER="${WORKER:-${SCRIPTS_DIR}/run_hyphy_relax_absrel.sh}"

ONLY="${ONLY:-both}"
export ONLY
[ -n "${HYPHY:-}" ] && export HYPHY

HYPHY_FA_DIR="${PHENOTYPE_DIR}/hyphy_fa"
[ -d "${HYPHY_FA_DIR}" ] || { echo "Missing ${HYPHY_FA_DIR}"; exit 1; }
[ -s "${RELAX_TREE}" ] || { echo "Missing RELAX tree: ${RELAX_TREE} -- run make_hyphy_trees.py for ${PHENOTYPE} first"; exit 1; }
[ -s "${ABSREL_TREE}" ] || { echo "Missing aBSREL tree: ${ABSREL_TREE} -- run make_hyphy_trees.py for ${PHENOTYPE} first"; exit 1; }
[ -x "${WORKER}" ] || { echo "Missing/not executable WORKER: ${WORKER}"; exit 1; }

mkdir -p "${OUTDIR}" "${TASKLIST_DIR}"

mapfile -t OGS < <(
  for f in "${HYPHY_FA_DIR}"/*.codon.fasta; do
    [ -e "${f}" ] || continue
    base="$(basename "${f}")"
    echo "${base%%.*}"
  done | sort
)
N=${#OGS[@]}
[ "${N}" -gt 0 ] || { echo "No *.codon.fasta found in ${HYPHY_FA_DIR}"; exit 1; }

TASK_ID="${SLURM_ARRAY_TASK_ID}"
STEP="${SLURM_ARRAY_TASK_COUNT:-40}"
START=$((TASK_ID-1))

TASK_OG_LIST="${TASKLIST_DIR}/task${TASK_ID}.txt"
: > "${TASK_OG_LIST}"
for ((i=START; i<N; i+=STEP)); do
  echo "${OGS[$i]}" >> "${TASK_OG_LIST}"
done
n_task=$(wc -l < "${TASK_OG_LIST}")

echo "PHENOTYPE=${PHENOTYPE}  ONLY=${ONLY}"
echo "PHENOTYPE_DIR=${PHENOTYPE_DIR}"
echo "OUTDIR=${OUTDIR}"
echo "Total OGs (genome-wide): ${N}"
echo "Task ${TASK_ID}: ${n_task} OGs -> ${TASK_OG_LIST}"

bash "${WORKER}" "${PHENOTYPE_DIR}" "${RELAX_TREE}" "${ABSREL_TREE}" "${OUTDIR}" "${TASK_OG_LIST}"
