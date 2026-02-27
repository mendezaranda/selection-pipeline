#!/bin/bash
#SBATCH --job-name=bpp_nt
#SBATCH -o /home/dmendez/outdir/bpp_nt_%A_%a.out
#SBATCH -e /home/dmendez/errdir/bpp_nt_%A_%a.err
#SBATCH -t 24:00:00
#SBATCH --gres=localtmp:120G
#SBATCH -N1
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-user=Daniel.MendezAranda@mdc-berlin.de
#SBATCH --array=1-10

set -eo pipefail

THREADS="${SLURM_CPUS_PER_TASK:-16}"

RESULTS_DIR="/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/03b_orthofinder/out/Results_Feb12"
SCO_DIR="${RESULTS_DIR}/Single_Copy_Orthologue_Sequences_renamed"

OUTDIR="/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/08_bpp"
ALL_CDS="${OUTDIR}/ALL.primary.cds.fa"

# CDS-only per-OG script (align CDS directly, optional nucleotide trimming).
# Default assumes you install og_to_codon2_phy.sh as /home/dmendez/scripts/og_to_nt_phy.sh.
# You can override at submit time: sbatch --export=OG_SCRIPT=/path/to/og_to_codon2_phy.sh build_bpp_nt_cds_only.sh
OG_SCRIPT_DEFAULT="/home/dmendez/scripts/og_to_nt_phy.sh"
OG_SCRIPT="/home/dmendez/scripts/og_to_codon2_phy.sh"

# Locus controls (override at submit time if needed)
MIN_TAXA="${MIN_TAXA:-10}"
DO_TRIMAL="${DO_TRIMAL:-1}"          # 1: trim nucleotide alignment (trimAl), 0: no trimming
FILTER_STOPS="${FILTER_STOPS:-1}"    # 1: drop sequences with internal stop codons (translation check)
FILTER_FRAME="${FILTER_FRAME:-1}"    # 1: drop sequences with length%3!=0

mkdir -p "${OUTDIR}"

[ -s "${ALL_CDS}" ] || { echo "Missing ALL_CDS: ${ALL_CDS}"; exit 1; }
[ -x "${OG_SCRIPT}" ] || { echo "Missing/Not executable OG_SCRIPT: ${OG_SCRIPT}"; exit 1; }

# OG list (protein fasta; used only to get member IDs per OG)
mapfile -t OGFILES < <(ls "${SCO_DIR}"/OG*.fa "${SCO_DIR}"/OG*.fasta 2>/dev/null | sort)
N=${#OGFILES[@]}
[ "${N}" -gt 0 ] || { echo "No OG*.fa found in ${SCO_DIR}"; exit 1; }

TASK_ID="${SLURM_ARRAY_TASK_ID}"

# Use Slurm array size as the stride (so changing --array=1-X doesn't require editing this script)
STEP="${SLURM_ARRAY_TASK_COUNT:-10}"
START=$((TASK_ID-1))

echo "Total OG files: ${N}"
echo "Task ${TASK_ID}: processing indices ${START}, ${START}+${STEP}, ..."

for ((i=START; i<N; i+=STEP)); do
  og="${OGFILES[$i]}"
  echo "Running OG: ${og}"

  MIN_TAXA="${MIN_TAXA}" DO_TRIMAL="${DO_TRIMAL}" FILTER_STOPS="${FILTER_STOPS}" FILTER_FRAME="${FILTER_FRAME}" \
    bash "${OG_SCRIPT}" "${og}" "${ALL_CDS}" "${OUTDIR}" "${THREADS}"
done
