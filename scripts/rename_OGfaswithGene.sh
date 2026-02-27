#!/bin/bash

OUTDIR=/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/09_pain_genes/out_2
MAP=/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/09_pain_genes/out_2/mouse_targets_to_OG.tsv

mkdir -p "$OUTDIR/og_fastas_by_gene"

python3 - <<'PY'
import csv, re
from pathlib import Path

outdir = Path("/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/09_pain_genes/out_2")
map_tsv = Path("/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/09_pain_genes/out_2/mouse_targets_to_OG.tsv")
src_dir = outdir / "og_fastas"
dst_dir = outdir / "og_fastas_by_gene"
dst_dir.mkdir(parents=True, exist_ok=True)

def clean(s):
    # safe filename
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', s)

made = 0
missing = 0

with map_tsv.open() as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for row in rd:
        if row.get("found","0") != "1":
            continue
        gene = row["mouse_gene"].strip() or "NA"
        og = row["orthogroup"].strip()
        if og == "NA" or not og:
            continue
        src = src_dir / f"{og}.fa"
        if not src.exists():
            missing += 1
            continue
        name = f"{clean(gene)}__{og}.fa"
        dst = dst_dir / name
        if dst.exists():
            continue
        # relative symlink for portability
        rel = Path("..") / "og_fastas" / src.name
        dst.symlink_to(rel)
        made += 1

print("symlinks made:", made)
print("missing OG fasta:", missing)
print("dst:", dst_dir)
PY
