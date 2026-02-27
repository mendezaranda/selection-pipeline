#!/usr/bin/env python3
# Compatible with Python 3.7–3.9

import argparse
import csv
import re
import sys
import shutil
from collections import defaultdict
from pathlib import Path

NPXP_RE = re.compile(r'([NX]P_\d+(?:\.\d+)?)')

def die(msg):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(2)

def sniff_delim(p):
    with Path(p).open("r", errors="ignore") as fh:
        line = fh.readline()
    if "\t" in line:
        return "\t"
    if "," in line:
        return ","
    return "\t"

def read_query_list(path):
    queries = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            queries.append(line)
    return queries

def normalize_accessions_from_token(token):
    out = set()
    m = NPXP_RE.search(token)
    if not m:
        return out
    acc = m.group(1)
    out.add(acc)
    out.add(acc.split(".")[0])
    return out

def load_protein2gene(tsv_path):
    delim = sniff_delim(tsv_path)

    with open(tsv_path, "r") as fh:
        rd = csv.DictReader(fh, delimiter=delim)
        if not rd.fieldnames:
            die("protein2gene file has no header: %s" % tsv_path)
        fields = [f.strip() for f in rd.fieldnames]

        prot_col = None
        gene_col = None
        for c in ["protein_id", "protein", "accession", "refseq", "Protein", "Accession"]:
            if c in fields:
                prot_col = c
                break
        for c in ["gene", "gene_name", "symbol", "Gene", "GeneName", "geneSymbol"]:
            if c in fields:
                gene_col = c
                break
        if prot_col is None:
            prot_col = fields[0]
        if gene_col is None:
            gene_col = fields[1] if len(fields) > 1 else fields[0]

    prot2gene = {}
    gene2prot = defaultdict(set)

    with open(tsv_path, "r") as fh:
        rd = csv.DictReader(fh, delimiter=delim)
        for row in rd:
            p = (row.get(prot_col) or "").strip()
            g = (row.get(gene_col) or "").strip()
            if not p or not g:
                continue
            m = NPXP_RE.search(p)
            if not m:
                continue
            acc = m.group(1)
            prot2gene[acc] = g
            prot2gene[acc.split(".")[0]] = g
            gene2prot[g].add(acc)
            gene2prot[g].add(acc.split(".")[0])

    return prot2gene, gene2prot, prot_col, gene_col

def parse_orthogroups_tsv(og_tsv, mouse_col):
    og2mouse_tokens = {}
    with open(og_tsv, "r") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if mouse_col not in header:
            die("mouse column '%s' not found. First header cols: %s"
                % (mouse_col, ", ".join(header[:20])))

        mi = header.index(mouse_col)

        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            og = parts[0].strip()
            if not og:
                continue
            cell = parts[mi].strip() if mi < len(parts) else ""
            if not cell:
                og2mouse_tokens[og] = []
                continue
            toks = [x.strip() for x in cell.split(",") if x.strip()]
            og2mouse_tokens[og] = toks

    return og2mouse_tokens

def build_mouse_accession_index(og2mouse_tokens):
    acc2ogs = defaultdict(set)
    for og, toks in og2mouse_tokens.items():
        for tok in toks:
            accs = normalize_accessions_from_token(tok)
            if accs:
                for a in accs:
                    acc2ogs[a].add(og)
            else:
                acc2ogs[tok].add(og)
    return acc2ogs

def find_og_fasta(og_fasta_dir, og):
    d = Path(og_fasta_dir)
    p1 = d / (og + ".fa")
    if p1.exists() and p1.stat().st_size > 0:
        return p1
    p2 = d / (og + ".fasta")
    if p2.exists() and p2.stat().st_size > 0:
        return p2
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orthogroups_tsv", required=True)
    ap.add_argument("--mouse_col", default="MusMus")
    ap.add_argument("--protein2gene_tsv", required=True)
    ap.add_argument("--query_list", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--og_fasta_dir", required=True,
                    help="Directory containing OG fasta files (e.g. .../Orthogroups/Orthogroup_Sequences or .../Single_Copy_Orthologue_Sequences_renamed)")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite copied fasta files if they exist")
    ap.add_argument("--out_og_list", default="", help="Optional: also write unique OG IDs here")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_fastas = outdir / "og_fastas"
    out_fastas.mkdir(parents=True, exist_ok=True)

    prot2gene, gene2prot, prot_col, gene_col = load_protein2gene(args.protein2gene_tsv)
    og2mouse_tokens = parse_orthogroups_tsv(args.orthogroups_tsv, args.mouse_col)
    acc2ogs = build_mouse_accession_index(og2mouse_tokens)
    queries = read_query_list(args.query_list)

    rows = []
    ogs = set()

    for q in queries:
        accs = normalize_accessions_from_token(q)
        if accs:
            acc = None
            for a in accs:
                if "." in a:
                    acc = a
                    break
            if acc is None:
                acc = sorted(accs)[0]

            gene = prot2gene.get(acc) or prot2gene.get(acc.split(".")[0]) or "NA"
            hits = set()
            for a in accs:
                hits |= acc2ogs.get(a, set())

            if not hits:
                rows.append((q, gene, "NA", "0"))
            else:
                for og in sorted(hits):
                    rows.append((q, gene, og, "1"))
                    ogs.add(og)
        else:
            gene = q
            prots = gene2prot.get(gene, set())
            if not prots:
                rows.append((q, gene, "NA", "0"))
                continue
            hits = set()
            for p in prots:
                hits |= acc2ogs.get(p, set())

            if not hits:
                rows.append((q, gene, "NA", "0"))
            else:
                for og in sorted(hits):
                    rows.append((q, gene, og, "1"))
                    ogs.add(og)

    # Write mapping TSV
    out_tsv = outdir / "mouse_targets_to_OG.tsv"
    with out_tsv.open("w") as oh:
        oh.write("query\tmouse_gene\torthogroup\tfound\n")
        for r in rows:
            oh.write("\t".join(r) + "\n")

    # Optional OG list
    if args.out_og_list:
        Path(args.out_og_list).write_text("\n".join(sorted(ogs)) + "\n")
    else:
        (outdir / "OGs_to_extract.txt").write_text("\n".join(sorted(ogs)) + "\n")

    # Copy OG fastas
    copied = 0
    missing = 0
    for og in sorted(ogs):
        src = find_og_fasta(args.og_fasta_dir, og)
        if src is None:
            missing += 1
            continue
        dst = out_fastas / (og + ".fa")
        if dst.exists() and (not args.overwrite):
            continue
        shutil.copyfile(str(src), str(dst))
        copied += 1

    sys.stderr.write("Wrote: %s\n" % out_tsv)
    sys.stderr.write("OGs: %d\n" % len(ogs))
    sys.stderr.write("Copied OG fastas: %d\n" % copied)
    sys.stderr.write("Missing OG fastas: %d (from %s)\n" % (missing, args.og_fasta_dir))
    sys.stderr.write("protein2gene columns used: %s -> %s\n" % (prot_col, gene_col))
    sys.stderr.write("OG fastas in: %s\n" % out_fastas)

if __name__ == "__main__":
    main()
