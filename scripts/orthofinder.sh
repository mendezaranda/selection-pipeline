#!/bin/bash
#
#SBATCH -n 32
#SBATCH -N 1
#SBATCH --mem-per-cpu=4G
#SBATCH -t 24:00:00
#SBATCH --job-name=orthofinder
#SBATCH -o /home/dmendez/outdir/%j_%x
#SBATCH -e /home/dmendez/errdir/%j_%x
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=Daniel.MendezAranda@mdc-berlin.de

# ================================================================
# Orthofinder 
# Daniel Mendez Aranda, 2025-11
# ================================================================

# --- USER SETTINGS -------------------------------------------------
source /home/dmendez/.bashrc
conda activate /fast/AG_Lewin/dmendez/.conda/envs/orthofinder
module load mpi

#cd /fast/AG_Lewin/dmendez/orthofinder/primary_transcripts
cd /fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/01b_proteomes

#for f in *fa ; do python /fast/AG_Lewin/dmendez/.conda/envs/orthofinder/bin/primary_transcript.py $f ; done

orthofinder -f /fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/01b_proteomes \
  -t 32 -a 32 \
  -M msa -S diamond  -o /fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/03b_orthofinder/out

conda deactivate 
date