#!/bin/bash
#SBATCH --job-name=absrel_fg
#SBATCH -o /home/dmendez/outdir/absrel_%A_%a.out
#SBATCH -e /home/dmendez/errdir/absrel_%A_%a.err
#SBATCH -t 24:00:00
#SBATCH -N1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --array=1-10%10

set -euo pipefail
source /home/dmendez/.bashrc

THREADS="${SLURM_CPUS_PER_TASK:-16}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

ALN_DIR="${ALN_DIR:-/fast/AG_Lewin/dmendez/selection_from_SCO/hyphy_fa}"
TREE_FG="${TREE_FG:-/fast/AG_Lewin/dmendez/species.hyphy_fg.nwk}"
OUTDIR="${OUTDIR:-/fast/AG_Lewin/dmendez/hyphy_absrel_fg}"
mkdir -p "${OUTDIR}"

command -v hyphy >/dev/null 2>&1 || { echo "hyphy not found in PATH"; exit 2; }
[ -s "${TREE_FG}" ] || { echo "Missing tree: ${TREE_FG}"; exit 2; }

mapfile -t ALNS < <(ls "${ALN_DIR}"/*.fasta "${ALN_DIR}"/*.fa 2>/dev/null | sort)
N=${#ALNS[@]}
[ "${N}" -gt 0 ] || { echo "No FASTA alignments found in ${ALN_DIR}"; exit 1; }

TASK_ID="${SLURM_ARRAY_TASK_ID}"
STEP="${SLURM_ARRAY_TASK_COUNT:-10}"
START=$((TASK_ID-1))

echo "Total alignments: ${N}"
echo "Task ${TASK_ID}: indices ${START}, ${START}+${STEP}, ..."

for ((i=START; i<N; i+=STEP)); do
  aln="${ALNS[$i]}"
  base="$(basename "$aln")"
  og="${base%%.*}"

  out_json="${OUTDIR}/${og}.absrel.json"
  out_txt="${OUTDIR}/${og}.absrel.log"

  if [ -s "${out_json}" ]; then
    echo "SKIP ${og} (exists)"
    continue
  fi

  echo "Running aBSREL: ${og}"

  work="${OUTDIR}/tmp_${SLURM_JOB_ID}_${og}"
  mkdir -p "${work}"
  cp -f "${aln}" "${work}/alignment.fasta"
  cp -f "${TREE_FG}" "${work}/tree_fg.nwk"

  (
    cd "${work}"
    hyphy absrel --alignment alignment.fasta --tree tree_fg.nwk --output result.json > run.log 2>&1
  ) || true

  if [ -s "${work}/result.json" ]; then
    mv -f "${work}/result.json" "${out_json}"
  else
    j="$(ls -1 "${work}"/*.json 2>/dev/null | head -n1 || true)"
    if [ -n "${j}" ] && [ -s "${j}" ]; then
      mv -f "${j}" "${out_json}"
    fi
  fi

  if [ -s "${work}/run.log" ]; then
    mv -f "${work}/run.log" "${out_txt}"
  fi

  rm -rf "${work}"

  if [ ! -s "${out_json}" ]; then
    echo "FAILED ${og} (no JSON produced). See ${out_txt}"
  else
    echo "OK ${og} -> ${out_json}"
  fi
done
