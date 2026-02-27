#!/usr/bin/env python3
import re
import sys
from collections import defaultdict

if len(sys.argv) != 3:
    sys.stderr.write("Usage: mikado_primary_gff.py <in.gff3> <out.primary.gff3>\n")
    sys.exit(2)

inp, outp = sys.argv[1], sys.argv[2]

id_re = re.compile(r'(?:^|;)ID=([^;]+)')
parent_re = re.compile(r'(?:^|;)Parent=([^;]+)')
# mikado mRNA IDs commonly like: <gene>.<iso>  (e.g., mikado.scaffold_1G4.1)
iso_re = re.compile(r'^(?P<gene>.+)\.(?P<iso>\d+)$')

# accept both mRNA and transcript as "transcript features"
TX_TYPES = {"mRNA", "transcript"}

# pass 1: collect transcripts by gene and detect which are explicitly primary
tx_by_gene = defaultdict(list)     # gene -> list of (iso, tx_id, is_primary)
all_tx = set()

def split_gff(line):
    # robust split: prefer tab, else any whitespace
    parts = line.rstrip("\n").split("\t")
    if len(parts) >= 9:
        return parts
    parts = re.split(r"\s+", line.rstrip("\n"), maxsplit=8)
    return parts if len(parts) >= 9 else None

with open(inp) as fh:
    for line in fh:
        if not line or line.startswith("#"):
            continue
        parts = split_gff(line)
        if not parts:
            continue
        ftype = parts[2]
        attrs = parts[8]
        if ftype not in TX_TYPES:
            continue

        m = id_re.search(attrs)
        if not m:
            continue
        tx_id = m.group(1)
        all_tx.add(tx_id)

        # infer gene id from Parent= if present, else from tx_id by iso pattern
        pm = parent_re.search(attrs)
        gene_id = pm.group(1).split(",")[0] if pm else None

        iso = None
        mm = iso_re.match(tx_id)
        if mm:
            iso = int(mm.group("iso"))
            gene_guess = mm.group("gene")
            if gene_id is None:
                gene_id = gene_guess

        is_primary = ("primary=True" in attrs)

        if gene_id is None:
            # if we cannot determine gene, treat tx as its own gene bucket
            gene_id = tx_id
            iso = iso if iso is not None else 10**9

        if iso is None:
            iso = 10**9  # unknown iso -> worst, unless primary=True

        tx_by_gene[gene_id].append((iso, tx_id, is_primary))

# choose one transcript per gene:
# 1) primary=True if any
# 2) else iso==1 if exists
# 3) else smallest iso
chosen_tx = set()
for gene, lst in tx_by_gene.items():
    prim = [x for x in lst if x[2]]
    if prim:
        # if multiple "primary=True", keep smallest iso among them
        pick = sorted(prim, key=lambda x: x[0])[0][1]
    else:
        iso1 = [x for x in lst if x[0] == 1]
        if iso1:
            pick = iso1[0][1]
        else:
            pick = sorted(lst, key=lambda x: x[0])[0][1]
    chosen_tx.add(pick)

# pass 2: write:
# - gene lines (optional: you can keep all, but this keeps only those that have at least one chosen transcript)
# - chosen transcript lines
# - child features whose Parent intersects chosen_tx
chosen_genes = set(tx_by_gene.keys())

def parents_in(attrs):
    m = parent_re.search(attrs)
    if not m:
        return []
    return m.group(1).split(",")

with open(inp) as fh, open(outp, "w") as oh:
    for line in fh:
        if line.startswith("#"):
            oh.write(line)
            continue
        parts = split_gff(line)
        if not parts:
            continue
        ftype = parts[2]
        attrs = parts[8]

        if ftype == "gene":
            # keep gene if it has any transcript bucket (conservative)
            gm = id_re.search(attrs)
            if gm and gm.group(1) in chosen_genes:
                oh.write(line)
            continue

        if ftype in TX_TYPES:
            m = id_re.search(attrs)
            if m and m.group(1) in chosen_tx:
                oh.write(line)
            continue

        # child features
        for p in parents_in(attrs):
            if p in chosen_tx:
                oh.write(line)
                break

sys.stderr.write(f"Genes with transcripts: {len(tx_by_gene)}\n")
sys.stderr.write(f"Chosen transcripts:     {len(chosen_tx)}\n")
sys.stderr.write(f"Wrote: {outp}\n")
