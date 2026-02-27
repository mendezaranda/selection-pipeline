#!/usr/bin/env python3

import argparse, re
from pathlib import Path
from collections import Counter

def read_ids_from_og_fa(p: Path):
    ids=[]
    with p.open() as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids

def fasta_headers_set(p: Path):
    s=set()
    with p.open() as fh:
        for line in fh:
            if line.startswith(">"):
                s.add(line[1:].split()[0])
    return s

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--sco_dir", required=True)
    ap.add_argument("--all_cds", required=True)
    ap.add_argument("--missing_list", required=True)
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--min_taxa", type=int, default=10)
    args=ap.parse_args()

    sco=Path(args.sco_dir)
    all_cds=Path(args.all_cds)
    missing=Path(args.missing_list)

    if not all_cds.exists(): raise SystemExit("missing ALL_CDS")
    all_ids = fasta_headers_set(all_cds)

    out=Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)

    counts=Counter()

    with out.open("w") as oh:
        oh.write("OG\tn_ids\tn_found\tstatus\tnote\n")
        for og in missing.read_text().splitlines():
            og=og.strip()
            if not og: continue
            ogfa = sco / f"{og}.fa"
            if not ogfa.exists():
                oh.write(f"{og}\t0\t0\tMISSING_OG_FASTA\t-\n")
                counts["MISSING_OG_FASTA"] += 1
                continue
            ids = read_ids_from_og_fa(ogfa)
            want=set(ids)
            found = sum(1 for x in want if x in all_ids)
            if found < args.min_taxa:
                oh.write(f"{og}\t{len(ids)}\t{found}\tTOO_FEW_FOUND\tMIN_TAXA={args.min_taxa}\n")
                counts["TOO_FEW_FOUND"] += 1
            else:
                oh.write(f"{og}\t{len(ids)}\t{found}\tFOUND_OK\t-\n")
                counts["FOUND_OK"] += 1

    print("Counts:")
    for k,v in counts.most_common():
        print(v, k)

if __name__=="__main__":
    main()
