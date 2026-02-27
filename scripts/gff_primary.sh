#!/bin/bash
gff="/fast/AG_Lewin/dmendez/genomes/BatSui/BatSui_Dustin/BatSui.edited.geneSymbol.gff"
primary_gff="/fast/AG_Lewin/dmendez/genomes/BatSui/BatSui_Dustin/BatSui.primary.gff"

python3 - <<'PY'
import re

inp  = "/fast/AG_Lewin/dmendez/genomes/BatSui/BatSui_Dustin/BatSui.edited.geneSymbol.gff"
outp = "/fast/AG_Lewin/dmendez/genomes/BatSui/BatSui_Dustin/BatSui.primary.gff"

primary_mrna = set()

# first pass: collect primary mRNA IDs
with open(inp) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        if parts[2] == "mRNA" and "primary=True" in parts[8]:
            m = re.search(r'ID=([^;]+)', parts[8])
            if m:
                primary_mrna.add(m.group(1))

# second pass: keep gene lines, primary mRNA lines, and features whose Parent is a primary mRNA
with open(inp) as fh, open(outp, "w") as oh:
    for line in fh:
        if line.startswith("#"):
            oh.write(line)
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        ftype = parts[2]
        attrs = parts[8]

        if ftype == "gene":
            oh.write(line)
            continue

        if ftype == "mRNA":
            m = re.search(r'ID=([^;]+)', attrs)
            if m and m.group(1) in primary_mrna:
                oh.write(line)
            continue

        # child features (CDS/exon/UTR etc)
        m = re.search(r'Parent=([^;]+)', attrs)
        if m:
            parent = m.group(1)
            if parent in primary_mrna:
                oh.write(line)

print(f"Primary mRNAs kept: {len(primary_mrna)}")
print(f"Wrote: {outp}")
PY
