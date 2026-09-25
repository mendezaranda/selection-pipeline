#!/usr/bin/env bash
# Copy this file to configs/env.sh (gitignored -- keep cluster-specific paths
# out of the repo) and edit the values for your run. Then, before invoking
# any script under scripts/paml/ or scripts/hyphy/:
#
#   source configs/env.sh
#
# Every script in scripts/paml/ and scripts/hyphy/ reads these as environment
# variables with a "${VAR:?...}" or "${VAR:-default}" pattern -- nothing
# parses this file directly, it's just a convenient place to set them all at
# once instead of exporting each one by hand. Every value can also be
# overridden per-invocation, e.g.:
#
#   RUN_BASE=/other/run/dir sbatch --export=ALL scripts/paml/build_codon_alignments_array.sh
#
# Scripts locate EACH OTHER (e.g. an array driver finding its per-OG worker
# script) relative to their own location on disk, not a hardcoded path --
# so scripts/paml/ and scripts/hyphy/ each work as a self-contained
# "scripts directory" wherever you clone or copy them to (a cluster's
# ~/scripts/, a flat directory you scp individual files into, etc.).

# ---- Run root ---------------------------------------------------------
# Everything else below is derived from this by default. This is the one
# variable you almost certainly need to change.
export RUN_BASE="/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01"

# ---- Standard subdirectories under RUN_BASE ----------------------------
# Override any single one of these if your layout differs from the default
# convention (RUN_BASE/<stage>/...).
export ORTHOFINDER_RESULTS_DIR="${ORTHOFINDER_RESULTS_DIR:-${RUN_BASE}/03b_orthofinder/out/Results_Feb12}"
export SCO_DIR="${SCO_DIR:-${ORTHOFINDER_RESULTS_DIR}/Single_Copy_Orthologue_Sequences_renamed}"
export ALL_CDS="${ALL_CDS:-${RUN_BASE}/08_bpp/ALL.primary.cds.fa}"
export PAML_DIR="${PAML_DIR:-${RUN_BASE}/06_paml}"
export HYPHY_DIR="${HYPHY_DIR:-${RUN_BASE}/07_hyphy}"

# ---- External tools -----------------------------------------------------
# Defaults assume they're on PATH (e.g. after `conda activate` /
# `module load`). Override with an absolute path otherwise.
export CODEML="${CODEML:-codeml}"
export HYPHY="${HYPHY:-hyphy}"
export TRIMAL="${TRIMAL:-trimal}"

# ---- Gene-symbol annotation ---------------------------------------------
# Genome-wide OG -> MusMus protein -> gene-symbol lookup, built ONCE (see
# scripts/annotation/extract_og_to_musmus_protein.sh) and reused by every
# phenotype and every test (PAML and HyPhy alike) -- not per-phenotype.
export OG2ACC_TSV="${OG2ACC_TSV:-${PAML_DIR}/OG_to_MusMus_protein.tsv}"
export ACC2GENE_TSV="${ACC2GENE_TSV:-${PAML_DIR}/MusMus_protein2gene.tsv}"

# ---- SLURM log directories -----------------------------------------------
# Referenced by the #SBATCH -o / -e headers in scripts/paml/*array*.sh and
# scripts/hyphy/*array*.sh. Those headers are static (SLURM reads them
# before any script code runs, so they can't reference a shell variable at
# submit time) -- create these two directories once, or edit the #SBATCH
# lines directly if you want your logs somewhere else.
mkdir -p "${HOME}/outdir" "${HOME}/errdir"
