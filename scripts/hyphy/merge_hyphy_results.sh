#!/usr/bin/env bash
#SBATCH --job-name=hyphy_merge
#SBATCH -o /home/%u/outdir/hyphy_merge_%j.out
#SBATCH -e /home/%u/errdir/hyphy_merge_%j.err
#SBATCH -t 00:30:00
#SBATCH -N1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -eo pipefail

# "Merge" step for the HyPhy RELAX/aBSREL genome-wide scan. Unlike the PAML
# array drivers (which concatenate per-task TSVs, unreliably, inside task 1
# of the array job itself -- see the caveat in build_codon_alignments_array.sh),
# HyPhy has nothing to concatenate: run_hyphy_relax_absrel.sh writes one
# JSON per OG directly into the shared OUTDIR, and parse_hyphy_results.py
# already aggregates across ALL of them by globbing that directory,
# regardless of which array task produced which file. So the correct merge
# is just: run parse_hyphy_results.py once, AFTER every array task has
# actually finished -- submit THIS as a separate job with a SLURM
# dependency on the array job, e.g.:
#
#   array_jid=$(sbatch --parsable --export=ALL,PHENOTYPE=only_HelKap \
#       run_hyphy_relax_absrel_array.sh)
#   sbatch --dependency=afterany:${array_jid} --export=ALL,PHENOTYPE=only_HelKap \
#       merge_hyphy_results.sh
#
# (afterany, not afterok -- a handful of individual OGs failing inside the
# array, which run_hyphy_relax_absrel.sh already tolerates and logs per-OG,
# shouldn't block building the summary/plots for everything that DID
# succeed. Check the array job's .err logs separately for per-OG failures.)
#
# Runs the full per-phenotype reporting pipeline in one go, mirroring the
# PAML pipeline (build_significance_table.py -> annotate_with_gene_names.py
# -> plot_paml_summary.py) with HyPhy's own parse + two-test-aware
# significance/plot scripts:
#   1. parse_hyphy_results.py            -- per-OG JSONs -> one combined TSV
#   2. build_hyphy_significance_table.py -- BH-FDR + Holm, RELAX and aBSREL
#                                            corrected independently
#   3. annotate_with_gene_names.py       -- adds gene_symbol (same script as
#                                            the PAML side -- genome-wide
#                                            OG->gene map, reused unchanged)
#   4. plot_hyphy_summary.py             -- run twice: --test relax, --test absrel
#
# source configs/env.sh (see configs/env.example.sh) first, or export
# RUN_BASE yourself, then set PHENOTYPE.

PHENOTYPE="${PHENOTYPE:?Set PHENOTYPE, e.g. only_HelKap | HelKap_HetGla | capsaicin | acid | acid_wBat}"
RUN_BASE="${RUN_BASE:?Set RUN_BASE (see configs/env.example.sh)}"
PAML_DIR="${PAML_DIR:-${RUN_BASE}/06_paml}"
HYPHY_DIR="${HYPHY_DIR:-${RUN_BASE}/07_hyphy}"

HYPHY_BASE="${HYPHY_BASE:-${HYPHY_DIR}/${PHENOTYPE}}"
OUTDIR="${OUTDIR:-${HYPHY_BASE}/hyphy_runs}"
ABSREL_TREE="${ABSREL_TREE:-${HYPHY_BASE}/hyphy_trees/${PHENOTYPE}.aBSREL_fg.nwk}"
PLOT_DIR="${PLOT_DIR:-${HYPHY_BASE}/plots}"

# Genome-wide gene-symbol lookup, built once by scripts/annotation/, reused
# across every phenotype/test (see configs/env.example.sh).
OG2ACC_TSV="${OG2ACC_TSV:-${PAML_DIR}/OG_to_MusMus_protein.tsv}"
ACC2GENE_TSV="${ACC2GENE_TSV:-${PAML_DIR}/MusMus_protein2gene.tsv}"

# Self-locating, like every other array/driver script in this repo -- finds
# its sibling scripts (parse_hyphy_results.py etc.) next to itself, and the
# PAML-side annotate_with_gene_names.py one directory over.
SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANNOTATION_DIR="${ANNOTATION_DIR:-${SCRIPTS_DIR}/../annotation}"

RESULTS_TSV="${OUTDIR}/${PHENOTYPE}_hyphy_relax_absrel_results.tsv"
SIG_TSV="${OUTDIR}/${PHENOTYPE}_hyphy_significance.tsv"
ANNOTATED_TSV="${OUTDIR}/${PHENOTYPE}_hyphy_significance.annotated.tsv"

[ -d "${OUTDIR}/relax" ] || [ -d "${OUTDIR}/absrel" ] || { echo "Neither ${OUTDIR}/relax nor ${OUTDIR}/absrel exists -- has run_hyphy_relax_absrel(_array).sh run for ${PHENOTYPE} yet?"; exit 1; }
[ -s "${ABSREL_TREE}" ] || { echo "Missing ${ABSREL_TREE}"; exit 1; }

echo "=== 1/4: parse_hyphy_results.py ==="
python3 "${SCRIPTS_DIR}/parse_hyphy_results.py" "${OUTDIR}" "${ABSREL_TREE}" "${RESULTS_TSV}"

echo
echo "=== 2/4: build_hyphy_significance_table.py ==="
python3 "${SCRIPTS_DIR}/build_hyphy_significance_table.py" "${RESULTS_TSV}" "${SIG_TSV}"

echo
echo "=== 3/4: annotate_with_gene_names.py ==="
if [ -s "${OG2ACC_TSV}" ] && [ -s "${ACC2GENE_TSV}" ]; then
  python3 "${ANNOTATION_DIR}/annotate_with_gene_names.py" "${SIG_TSV}" "${OG2ACC_TSV}" "${ACC2GENE_TSV}" "${ANNOTATED_TSV}"
else
  echo "SKIPPING annotation: OG2ACC_TSV (${OG2ACC_TSV}) or ACC2GENE_TSV (${ACC2GENE_TSV}) not found."
  echo "Set those env vars to the real paths (see configs/env.example.sh) and rerun, or annotate ${SIG_TSV} manually."
  cp "${SIG_TSV}" "${ANNOTATED_TSV}"
fi

echo
echo "=== 4/4: plot_hyphy_summary.py (relax + absrel, BH) ==="
mkdir -p "${PLOT_DIR}"
for test in relax absrel; do
  python3 "${SCRIPTS_DIR}/plot_hyphy_summary.py" "${ANNOTATED_TSV}" \
    --phenotype "${PHENOTYPE}" --test "${test}" --correction BH --outdir "${PLOT_DIR}" \
    || echo "  (skipped ${test} plot -- see message above, e.g. that test hasn't been run yet for this phenotype)"
done

echo
echo "Done. Outputs:"
echo "  ${RESULTS_TSV}"
echo "  ${SIG_TSV}"
echo "  ${ANNOTATED_TSV}"
echo "  ${PLOT_DIR}/${PHENOTYPE}_hyphy_{relax,absrel}_summary_BH.png"
