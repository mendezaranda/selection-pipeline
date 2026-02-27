#!/usr/bin/env python3
import re
import sys
from pathlib import Path

tx_re = re.compile(r"^(?P<gene>.+)\.t(?P<t>\d+)$")

def fasta_iter(p: Path):
    h=None; s=[]
    with p.open() as fh:
        for line in fh:
            line=line.rstrip("\n")
            if line.startswith(">"):
                if h: yield h,"".join(s)
                h=line; s=[]
            else:
                s.append(line.strip())
        if h: yield h,"".join(s)

def load_fasta_dict(p: Path):
    d={}
    for h,seq in fasta_iter(p):
        sid=h[1:].split()[0]
        d[sid]=seq
    return d

def choose_ids(ids):
    by_gene={}
    for sid in ids:
        m=tx_re.match(sid)
        if not m:
            # no .tN pattern; treat as its own gene (keep as-is)
            by_gene.setdefault(sid, []).append((0, sid))
            continue
        gene=m.group("gene")
        tn=int(m.group("t"))
        by_gene.setdefault(gene, []).append((tn, sid))

    chosen=set()
    for gene, lst in by_gene.items():
        # prefer t1 if present, else smallest tN
        lst_sorted=sorted(lst, key=lambda x: x[0])
        pick = None
        for tn,sid in lst_sorted:
            if tn == 1:
                pick = sid
                break
        if pick is None:
            pick = lst_sorted[0][1]
        chosen.add(pick)
    return chosen

def write_fasta(outp: Path, seqs: dict, keep_ids: set):
    with outp.open("w") as oh:
        for sid in sorted(keep_ids):
            seq = seqs[sid]
            oh.write(f">{sid}\n")
            for i in range(0,len(seq),60):
                oh.write(seq[i:i+60]+"\n")

def main():
    if len(sys.argv) != 5:
        print("Usage: keep_one_tx_per_gene.py <cds.fa> <pep.fa> <out.cds.fa> <out.pep.fa>", file=sys.stderr)
        sys.exit(2)

    cds_in = Path(sys.argv[1])
    pep_in = Path(sys.argv[2])
    cds_out = Path(sys.argv[3])
    pep_out = Path(sys.argv[4])

    cds = load_fasta_dict(cds_in)
    pep = load_fasta_dict(pep_in)

    chosen = choose_ids(set(cds.keys()) | set(pep.keys()))

    missing_in_cds = sorted([x for x in chosen if x not in cds])
    missing_in_pep = sorted([x for x in chosen if x not in pep])
    if missing_in_cds or missing_in_pep:
        print("ERROR: chosen IDs missing in one of the files.", file=sys.stderr)
        if missing_in_cds:
            print("Missing in CDS (first 10):", missing_in_cds[:10], file=sys.stderr)
        if missing_in_pep:
            print("Missing in PEP (first 10):", missing_in_pep[:10], file=sys.stderr)
        sys.exit(1)

    write_fasta(cds_out, cds, chosen)
    write_fasta(pep_out, pep, chosen)

    print("Input CDS:", len(cds), "Input PEP:", len(pep))
    print("Kept genes/transcripts:", len(chosen))
    print("Wrote:", cds_out)
    print("Wrote:", pep_out)

if __name__ == "__main__":
    main()
