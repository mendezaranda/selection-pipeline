#!/bin/bash
#SBATCH --job-name=bpp_nt
#SBATCH -o /home/dmendez/outdir/bpp_nt_%A_%a.out
#SBATCH -e /home/dmendez/errdir/bpp_nt_%A_%a.err
#SBATCH -t 24:00:00
#SBATCH --gres=localtmp:200G
#SBATCH -N1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --array=1-10

set -euo pipefail
THREADS="${SLURM_CPUS_PER_TASK:-8}"

RESULTS_DIR="/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/03b_orthofinder/out/Results_Feb12"
SCO_DIR="${RESULTS_DIR}/Single_Copy_Orthologue_Sequences_renamed"

OUTDIR="/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/08_bpp"
ALL_CDS="${OUTDIR}/ALL.primary.cds.fa"

OG_SCRIPT="/home/dmendez/scripts/og_to_codon2_phy.sh"

# >>> IMPORTANT TUNABLES <<<
# If you have e.g. 8 species total, set MIN_TAXA=8 (or lower).
MIN_TAXA="${MIN_TAXA:-10}"

# For transcriptome-heavy sets, try disabling filters first:
# FILTER_STOPS=0 FILTER_FRAME=0
DO_TRIMAL="${DO_TRIMAL:-1}"
FILTER_STOPS="${FILTER_STOPS:-0}"
FILTER_FRAME="${FILTER_FRAME:-0}"

mkdir -p "${OUTDIR}"
[ -s "${ALL_CDS}" ] || { echo "Missing ALL_CDS: ${ALL_CDS}"; exit 1; }
[ -x "${OG_SCRIPT}" ] || { echo "Missing/Not executable OG_SCRIPT: ${OG_SCRIPT}"; exit 1; }

mapfile -t OGFILES < <(ls "${SCO_DIR}"/OG*.fa "${SCO_DIR}"/OG*.fasta 2>/dev/null | sort)
N=${#OGFILES[@]}
[ "${N}" -gt 0 ] || { echo "No OG*.fa found in ${SCO_DIR}"; exit 1; }

TASK_ID="${SLURM_ARRAY_TASK_ID}"
STEP="${SLURM_ARRAY_TASK_COUNT:-10}"
START=$((TASK_ID-1))

echo "Total OG files: ${N}"
echo "Task ${TASK_ID}: processing indices ${START}, ${START}+${STEP}, ..."

for ((i=START; i<N; i+=STEP)); do
  og="${OGFILES[$i]}"
  MIN_TAXA="${MIN_TAXA}" DO_TRIMAL="${DO_TRIMAL}" FILTER_STOPS="${FILTER_STOPS}" FILTER_FRAME="${FILTER_FRAME}" \
    bash "${OG_SCRIPT}" "${og}" "${ALL_CDS}" "${OUTDIR}" "${THREADS}"
done

# Merge per-task status files into one (only from task 1 to avoid races)
if [ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]; then
  OUTDIR="/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/08_bpp"
  merged="${OUTDIR}/build_status.tsv"
  first=1
  rm -f "$merged"
  for f in "${OUTDIR}"/build_status.task*.tsv; do
    [ -s "$f" ] || continue
    if [ $first -eq 1 ]; then
      cat "$f" > "$merged"
      first=0
    else
      tail -n +2 "$f" >> "$merged"
    fi
  done
fi

echo "DONE task ${TASK_ID}"
