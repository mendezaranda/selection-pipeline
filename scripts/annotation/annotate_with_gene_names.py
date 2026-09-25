#!/usr/bin/env python3
"""
Adds a gene_symbol column to ANY OG-keyed table (a significance table, or the
raw codeml summary, or anything else with an "OG" column) by joining through
the mouse ortholog -- genome-wide, not just a curated candidate-gene list.

Chain: OG -> MusMus protein accession (OG_to_MusMus_protein.tsv) ->
       gene symbol (MusMus_protein2gene.tsv)

OG_to_MusMus_protein.tsv only needs to be generated ONCE, from
Single_Copy_Orthologue_Sequences_renamed, via extract_og_to_musmus_protein.sh
(already run) -- it covers every single-copy OG regardless of which
phenotype's codeml results you're annotating, so reuse the same file for all
of them; no need to regenerate per phenotype.

Usage:
    python3 annotate_with_gene_names.py <input.tsv> <OG_to_MusMus_protein.tsv> \
        <MusMus_protein2gene.tsv> <output_annotated.tsv>

<input.tsv>: any TSV with an "OG" column (e.g. a significance table from
    build_significance_table.py, for any phenotype).
"""
import sys
import csv


def main():
    if len(sys.argv) != 5:
        sys.exit(f"usage: {sys.argv[0]} <input.tsv> <OG_to_MusMus_protein.tsv> "
                  f"<MusMus_protein2gene.tsv> <output_annotated.tsv>")
    in_path, og2acc_path, acc2gene_path, out_path = sys.argv[1:5]

    og2acc = {}
    with open(og2acc_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            og2acc[row['OG']] = row['MusMus_protein_id']

    acc2gene = {}
    with open(acc2gene_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            acc2gene[row['protein_id']] = row['gene']

    og2gene = {og: acc2gene.get(acc, "NA") for og, acc in og2acc.items()}

    with open(in_path) as f:
        r = csv.DictReader(f, delimiter='\t')
        rows = list(r)
        fieldnames = r.fieldnames

    if "OG" not in fieldnames:
        sys.exit(f"'{in_path}' has no OG column (found: {fieldnames})")

    new_fieldnames = fieldnames[:1] + ["gene_symbol"] + fieldnames[1:]
    n_named = 0
    for row in rows:
        row["gene_symbol"] = og2gene.get(row["OG"], "NA")
        if row["gene_symbol"] != "NA":
            n_named += 1

    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=new_fieldnames, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)

    print(f"Rows: {len(rows)}, resolved a gene symbol for {n_named} ({100*n_named/len(rows):.1f}%)")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
