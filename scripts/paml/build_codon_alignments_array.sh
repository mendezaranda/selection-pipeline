#!/usr/bin/env bash
#SBATCH --job-name=paml_hyphy
#SBATCH -o /home/%u/outdir/paml_hyphy_%A_%a.out
#SBATCH -e /home/%u/errdir/paml_hyphy_%A_%a.err
#SBATCH -t 24:00:00
#SBATCH --gres=localtmp:200G
#SBATCH -N1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --array=1-20

set -eo pipefail
# source configs/env.sh (see configs/env.example.sh) before submitting, or
# export RUN_BASE / SCO_DIR / ALL_CDS / CODEML yourself.
#
# NOTE on the #SBATCH -o/-e paths above: SLURM reads those before any shell
# code runs, so they can't reference $HOME or another variable at submit
# time -- edit them directly if you want logs somewhere other than
# ~/outdir, ~/errdir (see configs/env.example.sh, which creates those two
# directories for you).

THREADS="${SLURM_CPUS_PER_TASK:-16}"

RUN_BASE="${RUN_BASE:?Set RUN_BASE (see configs/env.example.sh), e.g. export RUN_BASE=/fast/.../2026-02-10_run01}"
SCO_DIR="${SCO_DIR:-${RUN_BASE}/03b_orthofinder/out/Results_Feb12/Single_Copy_Orthologue_Sequences_renamed}"
ALL_CDS="${ALL_CDS:-${RUN_BASE}/08_bpp/ALL.primary.cds.fa}"
PAML_DIR="${PAML_DIR:-${RUN_BASE}/06_paml}"

# This is the FULL, from-scratch build for the baseline phenotype -- run it
# once (typically for your single simplest foreground, e.g. one focal
# species), then reuse its alignments for every other phenotype via
# run_codeml_from_existing_array.sh instead of rerunning mafft/trimAl from
# scratch each time.
PHENOTYPE="${PHENOTYPE:-only_HelKap}"
OUTDIR="${OUTDIR:-${PAML_DIR}/${PHENOTYPE}}"

# Foreground-labeled species tree for PAML (leaves to test marked #1)
TREE_FG="${TREE_FG:-${OUTDIR}/species.paml_fg.nwk}"

# Scripts locate each other relative to THIS script's own location, so
# scripts/paml/ works as a self-contained directory wherever you clone or
# copy it to (no hardcoded /home/<user>/scripts/ path).
SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OG_SCRIPT="${OG_SCRIPT:-${SCRIPTS_DIR}/og_to_codon_paml_codeml.sh}"

MIN_TAXA="${MIN_TAXA:-8}"
DO_TRIMAL="${DO_TRIMAL:-1}"
FILTER_STOPS="${FILTER_STOPS:-1}"
FILTER_FRAME="${FILTER_FRAME:-1}"
RUN_CODEML="${RUN_CODEML:-1}"
CODEML="${CODEML:-codeml}"

mkdir -p "${OUTDIR}"
[ -s "${ALL_CDS}" ] || { echo "Missing ALL_CDS: ${ALL_CDS}"; exit 1; }
[ -s "${TREE_FG}" ] || { echo "Missing TREE_FG: ${TREE_FG}"; exit 1; }
[ -x "${OG_SCRIPT}" ] || { echo "Missing/Not executable OG_SCRIPT: ${OG_SCRIPT}"; exit 1; }

mapfile -t OGFILES < <(ls "${SCO_DIR}"/OG*.fa "${SCO_DIR}"/OG*.fasta 2>/dev/null | sort)
N=${#OGFILES[@]}
[ "${N}" -gt 0 ] || { echo "No OG*.fa found in ${SCO_DIR}"; exit 1; }

TASK_ID="${SLURM_ARRAY_TASK_ID}"
STEP="${SLURM_ARRAY_TASK_COUNT:-40}"
START=$((TASK_ID-1))

echo "Total OG files: ${N}"
echo "Task ${TASK_ID}: processing indices ${START}, ${START}+${STEP}, ..."

for ((i=START; i<N; i+=STEP)); do
  og="${OGFILES[$i]}"
  MIN_TAXA="${MIN_TAXA}" DO_TRIMAL="${DO_TRIMAL}" FILTER_STOPS="${FILTER_STOPS}" FILTER_FRAME="${FILTER_FRAME}" RUN_CODEML="${RUN_CODEML}" CODEML="${CODEML}" \
    bash "${OG_SCRIPT}" "${og}" "${ALL_CDS}" "${OUTDIR}" "${TREE_FG}" "${THREADS}"
done

# Merge per-task summaries (do this only in task 1 to avoid races writing
# the SAME file at once -- this does NOT wait for other tasks to finish; if
# you need a guaranteed-complete merge, run this merge step as a separate
# job with --dependency=afterany:<this array job's id> instead, the same
# way scripts/hyphy/merge_hyphy_results.sh does it for the HyPhy side).
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
