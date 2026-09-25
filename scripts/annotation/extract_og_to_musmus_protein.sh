#!/usr/bin/env bash
# Extracts the MusMus NCBI protein accession from every single-copy orthogroup
# fasta, producing a 2-column OG_ID -> MusMus_protein_accession table.
# This can be joined locally against MusMus_protein2gene.tsv to get a gene
# name for every OG (not just the curated pain-gene candidates).
#
# Run this ONCE per OrthoFinder run -- the result is genome-wide and reused
# by every phenotype and every test (PAML and HyPhy alike), not regenerated
# per phenotype. See OG2ACC_TSV / SCO_DIR in configs/env.example.sh.
#
# Usage:
#   cd "${SCO_DIR}"   # Single_Copy_Orthologue_Sequences_renamed
#   bash extract_og_to_musmus_protein.sh > "${OG2ACC_TSV}"

set -u
echo -e "OG\tMusMus_protein_id"
for f in OG*.fa; do
  og="${f%.fa}"
  acc=$(grep -m1 "^>MusMus_" "$f" | sed 's/^>MusMus_//')
  if [ -z "${acc}" ]; then
    acc="NA"
  fi
  echo -e "${og}\t${acc}"
done
