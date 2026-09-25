#!/bin/bash

set -eo pipefail
source /home/dmendez/.bashrc

OG_FA="${1:?Usage: bash og_to_codon_paml_hyphy.sh <OGxxxx.fa> <ALL_CDS.fa> <OUTDIR> [THREADS]}"
ALL_CDS="${2:?Usage: bash og_to_codon_paml_hyphy.sh <OGxxxx.fa> <ALL_CDS.fa> <OUTDIR> [THREADS]}"
OUTDIR="${3:?Usage: bash og_to_codon_paml_hyphy.sh <OGxxxx.fa> <ALL_CDS.fa> <OUTDIR> [THREADS]}"
THREADS="${4:-8}"

MIN_TAXA="${MIN_TAXA:-10}"
DO_TRIMAL="${DO_TRIMAL:-1}"           # trim AA alignment before backtranslation
TRIMAL="${TRIMAL:-/fast/AG_Lewin/dmendez/.conda/envs/funannotate/bin/trimal}"
FILTER_STOPS="${FILTER_STOPS:-1}"    # drop internal stops in CDS (translation check)
FILTER_FRAME="${FILTER_FRAME:-1}"    # drop len%3!=0

command -v mafft >/dev/null 2>&1 || { echo "mafft not found"; exit 2; }
if [ "${DO_TRIMAL}" -eq 1 ]; then
  command -v "${TRIMAL}" >/dev/null 2>&1 || { echo "trimAl not found at ${TRIMAL}"; exit 2; }
fi

mkdir -p "${OUTDIR}/tmp" "${OUTDIR}/aa_aln" "${OUTDIR}/codon_aln" "${OUTDIR}/paml_phy" "${OUTDIR}/hyphy_fa" "${OUTDIR}/logs"

ogbase="$(basename "${OG_FA}")"
ogname="${ogbase%%.*}"
work="${OUTDIR}/tmp/${ogname}"
mkdir -p "${work}"

aa_clean="${work}/${ogname}.aa.clean.fa"
ids="${work}/${ogname}.ids"
cds_raw="${work}/${ogname}.cds.raw.fa"
cds_filt="${work}/${ogname}.cds.filt.fa"
aa_from_cds="${work}/${ogname}.aa.fromcds.fa"

aa_aln="${OUTDIR}/aa_aln/${ogname}.aa.aln.fa"
aa_trim="${OUTDIR}/aa_aln/${ogname}.aa.aln.trim.fa"
codon_fa="${OUTDIR}/codon_aln/${ogname}.codon.aln.fa"

paml_phy="${OUTDIR}/paml_phy/${ogname}.codon.phy"
hyphy_fa="${OUTDIR}/hyphy_fa/${ogname}.codon.fasta"
log="${OUTDIR}/logs/${ogname}.log"

# Skip if already done
if [ -s "${paml_phy}" ] && [ -s "${hyphy_fa}" ]; then
  echo "SKIP ${ogname} (exists)"
  exit 0
fi

# 1) Clean OG protein fasta headers to first token only; collect IDs
awk '
  BEGIN{FS=" "}
  /^>/{h=$1; sub(/^>/,"",h); print ">" h; next}
  {gsub(/[ \t\r]/,""); print}
' "${OG_FA}" > "${aa_clean}"
grep '^>' "${aa_clean}" | sed 's/^>//' > "${ids}"

# 2) Extract matching CDS from ALL_CDS
python3 - <<PY > "${log}" 2>&1
ids_file = "${ids}"
cds_all  = "${ALL_CDS}"
out_cds  = "${cds_raw}"
og = "${ogname}"

want = [x.strip() for x in open(ids_file) if x.strip()]
want_set = set(want)

def fasta_iter(fn):
    name=None; seq=[]
    for line in open(fn):
        line=line.rstrip("\n")
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name=line[1:].split()[0]
            seq=[]
        else:
            seq.append(line.strip())
    if name:
        yield name, "".join(seq)

found=[]
with open(out_cds,"w") as o:
    for name, seq in fasta_iter(cds_all):
        if name in want_set:
            seq = seq.upper().replace("U","T")
            o.write(f">{name}\n{seq}\n")
            found.append(name)

missing=[x for x in want if x not in set(found)]
print(f"{og}: requested {len(want)} IDs; found {len(found)} CDS; missing {len(missing)}")
if missing:
    print("Missing (first 10): " + ",".join(missing[:10]))
PY

n_cds="$(grep -c '^>' "${cds_raw}" || true)"
if [ "${n_cds}" -lt "${MIN_TAXA}" ]; then
  echo "SKIP ${ogname} (only ${n_cds} CDS found; MIN_TAXA=${MIN_TAXA})"
  exit 0
fi

# 3) Filter CDS (optional): len%3 and internal stops
python3 - <<PY >> "${log}" 2>&1
FILTER_STOPS = int("${FILTER_STOPS}")
FILTER_FRAME = int("${FILTER_FRAME}")

inp = "${cds_raw}"
out = "${cds_filt}"

genetic_code = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def read_fasta(fn):
    name=None; seq=[]
    for line in open(fn):
        line=line.strip()
        if not line: continue
        if line.startswith(">"):
            if name: yield name, "".join(seq)
            name=line[1:].split()[0]
            seq=[]
        else:
            seq.append(line)
    if name: yield name, "".join(seq)

def has_internal_stop(seq):
    codons=[seq[i:i+3] for i in range(0,len(seq)-2,3)]
    aas=[]
    for c in codons:
        if any(x not in "ACGT" for x in c): aas.append("X")
        else: aas.append(genetic_code.get(c,"X"))
    return "*" in aas[:-1]

kept=0
dropped=[]
with open(out,"w") as o:
    for name, seq in read_fasta(inp):
        seq=seq.upper().replace("U","T")
        if FILTER_FRAME and (len(seq)%3!=0):
            dropped.append((name,"len_not_mod3"))
            continue
        if FILTER_STOPS and has_internal_stop(seq):
            dropped.append((name,"internal_stop"))
            continue
        o.write(f">{name}\n{seq}\n")
        kept+=1

print(f"Filter kept {kept}, dropped {len(dropped)}")
if dropped:
    print("Dropped (first 10): " + ",".join([f"{n}:{r}" for n,r in dropped[:10]]))
PY

n_filt="$(grep -c '^>' "${cds_filt}" || true)"
if [ "${n_filt}" -lt "${MIN_TAXA}" ]; then
  echo "SKIP ${ogname} (after filtering only ${n_filt}; MIN_TAXA=${MIN_TAXA})"
  exit 0
fi

# 4) Translate filtered CDS to AA (for alignment guidance)
python3 - <<PY >> "${log}" 2>&1
inp="${cds_filt}"
out="${aa_from_cds}"

genetic_code = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def read_fasta(fn):
    name=None; seq=[]
    for line in open(fn):
        line=line.strip()
        if not line: continue
        if line.startswith(">"):
            if name: yield name, "".join(seq)
            name=line[1:].split()[0]
            seq=[]
        else:
            seq.append(line)
    if name: yield name, "".join(seq)

def translate(seq):
    aas=[]
    for i in range(0,len(seq)-2,3):
        c=seq[i:i+3]
        if any(x not in "ACGT" for x in c): aas.append("X")
        else: aas.append(genetic_code.get(c,"X"))
    # drop terminal stop if present
    if aas and aas[-1] == "*":
        aas = aas[:-1]
    return "".join(aas)

with open(out,"w") as o:
    for name, seq in read_fasta(inp):
        seq=seq.upper().replace("U","T")
        o.write(f">{name}\n{translate(seq)}\n")
PY

# 5) Align AA, optional trim
mafft --auto --thread "${THREADS}" "${aa_from_cds}" > "${aa_aln}"

if [ "${DO_TRIMAL}" -eq 1 ]; then
  "${TRIMAL}" -in "${aa_aln}" -out "${aa_trim}" -automated1
else
  cp "${aa_aln}" "${aa_trim}"
fi

# 6) Backtranslate AA alignment -> codon alignment (triplet gaps), then rename to species-only
python3 - <<PY >> "${log}" 2>&1
import sys

AA_ALN="${aa_trim}"
CDS="${cds_filt}"
OUT="${codon_fa}"
MIN_TAXA=int("${MIN_TAXA}")

genetic_code = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def read_fasta(fn):
    name=None; seq=[]
    for line in open(fn):
        line=line.strip()
        if not line: continue
        if line.startswith(">"):
            if name: yield name, "".join(seq)
            name=line[1:].split()[0]
            seq=[]
        else:
            seq.append(line)
    if name: yield name, "".join(seq)

def aa_from_codon(c):
    if any(x not in "ACGT" for x in c): return "X"
    return genetic_code.get(c,"X")

aa = {n:s for n,s in read_fasta(AA_ALN)}
cds = {n:s.upper().replace("U","T") for n,s in read_fasta(CDS)}

out={}
bad=[]
for n, aln in aa.items():
    if n not in cds:
        bad.append((n,"no_cds"))
        continue
    cds_seq = cds[n]
    aln_nogap = aln.replace("-","").replace(".","")
    need = len(aln_nogap)
    if len(cds_seq) < need*3:
        bad.append((n,"cds_too_short"))
        continue

    cds_trim = cds_seq[:need*3]
    codons = [cds_trim[i:i+3] for i in range(0,len(cds_trim),3)]

    # reject internal stops
    aa_check = "".join(aa_from_codon(c) for c in codons)
    if "*" in aa_check[:-1]:
        bad.append((n,"internal_stop"))
        continue

    cod_i=0
    cod_aln=[]
    mismatch=0
    for ch in aln:
        if ch in "-.":
            cod_aln.append("---")
        else:
            c = codons[cod_i]; cod_i += 1
            a = aa_from_codon(c)
            if ch != "X" and a != "X" and ch != a:
                mismatch += 1
            cod_aln.append(c)

    if mismatch:
        bad.append((n,f"aa_mismatch:{mismatch}"))
        continue

    # rename to species-only (assumes IDs are Species_something)
    sp = n.split("_",1)[0]
    out[sp] = "".join(cod_aln)

if len(out) < MIN_TAXA:
    open(OUT,"w").close()
    sys.exit(0)

L=None
for sp, seq in out.items():
    if L is None: L=len(seq)
    if len(seq)!=L:
        open(OUT,"w").close()
        sys.exit(0)

if L % 3 != 0:
    open(OUT,"w").close()
    sys.exit(0)

with open(OUT,"w") as o:
    for sp in sorted(out):
        o.write(f">{sp}\n{out[sp]}\n")

if bad:
    sys.stderr.write("Dropped (first 10): " + ",".join([f"{x}:{y}" for x,y in bad[:10]]) + "\n")
PY

if [ ! -s "${codon_fa}" ]; then
  echo "SKIP ${ogname} (no codon alignment produced; see ${log})"
  exit 0
fi

# 7) Write HyPhy FASTA (already FASTA) and PAML PHYLIP
cp "${codon_fa}" "${hyphy_fa}"

python3 - <<PY >> "${log}" 2>&1
inp="${codon_fa}"
out="${paml_phy}"

def read_fasta(fn):
    name=None; seq=[]
    for line in open(fn):
        line=line.strip()
        if not line: continue
        if line.startswith(">"):
            if name: yield name, "".join(seq)
            name=line[1:].split()[0]
            seq=[]
        else:
            seq.append(line)
    if name: yield name, "".join(seq)

recs=list(read_fasta(inp))
n=len(recs)
L=len(recs[0][1])
for nme, seq in recs:
    if len(seq)!=L:
        raise SystemExit("Length mismatch")

with open(out,"w") as o:
    o.write(f"{n} {L}\n")
    for nme, seq in recs:
        o.write(f"{nme}  {seq}\n")
PY

echo "OK ${ogname} -> ${paml_phy} and ${hyphy_fa}"
