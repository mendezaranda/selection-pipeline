#!/bin/bash
#SBATCH --job-name=bpp_missing20
#SBATCH -o /home/dmendez/outdir/bpp_missing20_%A_%a.out
#SBATCH -e /home/dmendez/errdir/bpp_missing20_%A_%a.err
#SBATCH -t 24:00:00
#SBATCH -N1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --gres=localtmp:100G
#SBATCH --array=1-20

set -euo pipefail

OUTDIR=/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/08_bpp
RESULTS_DIR=/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/03b_orthofinder/out/Results_Feb12
SCO_DIR=$RESULTS_DIR/Single_Copy_Orthologue_Sequences_renamed
ALL_CDS=$OUTDIR/ALL.primary.cds.fa
OG_SCRIPT=/home/dmendez/scripts/og_to_codon2_phy.sh
THREADS="${SLURM_CPUS_PER_TASK:-16}"

MISSING_LIST="$OUTDIR/missing_ogs.txt"

[ -s "$MISSING_LIST" ] || { echo "Missing list not found/empty: $MISSING_LIST"; exit 1; }
[ -s "$ALL_CDS" ] || { echo "Missing ALL_CDS: $ALL_CDS"; exit 1; }
[ -x "$OG_SCRIPT" ] || { echo "OG_SCRIPT not executable: $OG_SCRIPT"; exit 1; }
[ -d "$SCO_DIR" ] || { echo "SCO_DIR missing: $SCO_DIR"; exit 1; }

mapfile -t MOGS < "$MISSING_LIST"
N=${#MOGS[@]}

TASK_ID="${SLURM_ARRAY_TASK_ID}"
NTASKS="${SLURM_ARRAY_TASK_COUNT:-20}"
START=$((TASK_ID-1))
STEP=$NTASKS

echo "Missing OGs total: $N"
echo "Task $TASK_ID: indices $START, $START+$STEP, ..."

FAILLOG="${OUTDIR}/failed_missing.task${TASK_ID}.txt"
: > "$FAILLOG"

for ((i=START; i<N; i+=STEP)); do
  og="${MOGS[$i]}"
  ogfa="${SCO_DIR}/${og}.fa"
  if [ ! -s "$ogfa" ]; then
    echo -e "NOFA\t$og" >> "$FAILLOG"
    echo "Missing OG fasta: $ogfa"
    continue
  fi

  echo "Running OG: $og"
  bash "$OG_SCRIPT" "$ogfa" "$ALL_CDS" "$OUTDIR" "$THREADS" \
    || { rc=$?; echo -e "${rc}\t${og}" >> "$FAILLOG"; echo "FAILED rc=$rc: $og"; continue; }
done

echo "Done task $TASK_ID. Failures logged to $FAILLOG"
