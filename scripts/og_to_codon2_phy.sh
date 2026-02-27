#!/bin/bash

set -eo pipefail
source /home/dmendez/.bashrc
conda activate orthofinder

# Usage:
#   bash og_to_codon2_phy.sh <OGxxxx.fa> <ALL_CDS.fa> <OUTDIR> [THREADS]
OG_FA="${1:?Usage: bash og_to_codon2_phy.sh <OGxxxx.fa> <ALL_CDS.fa> <OUTDIR> [THREADS]}"
ALL_CDS="${2:?Usage: bash og_to_codon2_phy.sh <OGxxxx.fa> <ALL_CDS.fa> <OUTDIR> [THREADS]}"
OUTDIR="${3:?Usage: bash og_to_codon2_phy.sh <OGxxxx.fa> <ALL_CDS.fa> <OUTDIR> [THREADS]}"
THREADS="${4:-8}"

# Controls (override via env)
MIN_TAXA="${MIN_TAXA:-10}"            # IMPORTANT: set this appropriately
DO_TRIMAL="${DO_TRIMAL:-1}"
FILTER_STOPS="${FILTER_STOPS:-1}"
FILTER_FRAME="${FILTER_FRAME:-1}"

# Tools (conda-friendly)
MAFFT_BIN="${MAFFT_BIN:-mafft}"
TRIMAL_BIN="${TRIMAL_BIN:-trimal}"    # expect in PATH if DO_TRIMAL=1

#command -v "${MAFFT_BIN}" >/dev/null 2>&1 || { echo "mafft not found: ${MAFFT_BIN}"; exit 2; }
#if [ "${DO_TRIMAL}" -eq 1 ]; then
#  command -v "${TRIMAL_BIN}" >/dev/null 2>&1 || { echo "trimAl not found: ${TRIMAL_BIN}"; exit 2; }
#fi

mkdir -p "${OUTDIR}/tmp" "${OUTDIR}/nt_aln" "${OUTDIR}/phy_blocks" "${OUTDIR}/logs"

ogbase="$(basename "${OG_FA}")"
ogname="${ogbase%%.*}"
work="${OUTDIR}/tmp/${ogname}"
mkdir -p "${work}"

aa_clean="${work}/${ogname}.aa.fa"
ids="${work}/${ogname}.ids"
cds_raw="${work}/${ogname}.cds.raw.fa"
cds_filt="${work}/${ogname}.cds.filt.fa"

nt_aln="${OUTDIR}/nt_aln/${ogname}.nt.aln.fa"
nt_trim="${OUTDIR}/nt_aln/${ogname}.nt.aln.trim.fa"
phy_out="${OUTDIR}/phy_blocks/${ogname}.phy"
log="${OUTDIR}/logs/${ogname}.log"

# status file: one per array task (no locks)
TASK_TAG="${SLURM_ARRAY_TASK_ID:-0}"
STATUS_TSV="${OUTDIR}/build_status.task${TASK_TAG}.tsv"

append_status () {
  local og="$1" st="$2" nids="$3" nfound="$4" naf="$5" note="$6"
  if [ ! -s "${STATUS_TSV}" ]; then
    echo -e "OG\tstatus\tn_ids\tn_found\tn_after_filter\tnote" > "${STATUS_TSV}"
  fi
  echo -e "${og}\t${st}\t${nids}\t${nfound}\t${naf}\t${note}" >> "${STATUS_TSV}"
}

if [ -s "${phy_out}" ]; then
  append_status "${ogname}" "SKIP_EXISTS" "NA" "NA" "NA" "phy_exists"
  exit 0
fi

# 1) Clean OG protein fasta headers to first token only; collect IDs
awk '
  BEGIN{FS=" "}
  /^>/{h=$1; sub(/^>/, "", h); print ">" h; next}
  {gsub(/[ \t\r]/,""); print}
' "${OG_FA}" > "${aa_clean}"
grep '^>' "${aa_clean}" | sed 's/^>//' > "${ids}"
n_ids="$(wc -l < "${ids}" | tr -d ' ')"

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

found = []
with open(out_cds, "w") as o:
    for name, seq in fasta_iter(cds_all):
        if name in want_set:
            seq = seq.upper().replace("U","T")
            o.write(f">{name}\n{seq}\n")
            found.append(name)

missing = [x for x in want if x not in set(found)]
print(f"{og}: requested {len(want)} IDs; found {len(found)} CDS; missing {len(missing)}")
if missing:
    print("Missing (first 20): " + ",".join(missing[:20]))
PY

n_found="$(grep -c '^>' "${cds_raw}" || true)"
if [ "${n_found}" -lt "${MIN_TAXA}" ]; then
  append_status "${ogname}" "SKIP_TOO_FEW_FOUND" "${n_ids}" "${n_found}" "NA" "MIN_TAXA=${MIN_TAXA}"
  exit 0
fi

# 3) Optional CDS filtering (frame + internal stops)
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
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
    aas = []
    for c in codons:
        if any(x not in "ACGT" for x in c):
            aas.append("X")
        else:
            aas.append(genetic_code.get(c, "X"))
    return "*" in aas[:-1]   # terminal stop allowed

kept = 0
dropped = []
with open(out, "w") as o:
    for name, seq in read_fasta(inp):
        seq = seq.upper().replace("U","T")
        if FILTER_FRAME and (len(seq) % 3 != 0):
            dropped.append((name, "len_not_mod3"))
            continue
        if FILTER_STOPS and has_internal_stop(seq):
            dropped.append((name, "internal_stop"))
            continue
        o.write(f">{name}\n{seq}\n")
        kept += 1

print(f"Filter kept {kept}, dropped {len(dropped)}")
if dropped:
    print("Dropped (first 20): " + ",".join([f"{n}:{r}" for n,r in dropped[:20]]))
PY

n_after="$(grep -c '^>' "${cds_filt}" || true)"
if [ "${n_after}" -lt "${MIN_TAXA}" ]; then
  append_status "${ogname}" "SKIP_TOO_FEW_AFTER_FILTER" "${n_ids}" "${n_found}" "${n_after}" "FILTER_FRAME=${FILTER_FRAME};FILTER_STOPS=${FILTER_STOPS}"
  exit 0
fi

# 4) Align nucleotide CDS
"${MAFFT_BIN}" --auto --thread "${THREADS}" "${cds_filt}" > "${nt_aln}"

# 5) Trim (optional)
if [ "${DO_TRIMAL}" -eq 1 ]; then
  "${TRIMAL_BIN}" -in "${nt_aln}" -out "${nt_trim}" -automated1
	# after trimAl produced $nt_trim
	if [ ! -s "$nt_trim" ] || [ "$(grep -c '^>' "$nt_trim")" -lt "$MIN_TAXA" ]; then
		cp "$nt_aln" "$nt_trim"
	fi
else
  cp "${nt_aln}" "${nt_trim}"
fi

# 6) FASTA -> PHYLIP
python3 - <<PY >> "${log}" 2>&1
import sys
fa = "${nt_trim}"
out = "${phy_out}"
min_taxa = int("${MIN_TAXA}")

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

records=[]
for name, seq in read_fasta(fa):
    seq = seq.upper().replace("U","T")
    # species name is prefix before first underscore, if present; otherwise keep whole token
    sp = name.split("_", 1)[0]
    tip = f"{sp}^{sp}"
    records.append((tip, seq))

if len(records) < min_taxa:
    open(out, "w").close()
    sys.exit(0)

L = len(records[0][1])
if L == 0:
    open(out, "w").close()
    sys.exit(0)

for t,s in records:
    if len(s) != L:
        open(out, "w").close()
        sys.exit(0)

with open(out, "w") as o:
    o.write(f"{len(records)} {L}\n")
    for tip, seq in records:
        o.write(f"{tip}  {seq}\n")
PY

if [ ! -s "${phy_out}" ]; then
  append_status "${ogname}" "SKIP_NO_PHYLIP" "${n_ids}" "${n_found}" "${n_after}" "see_log"
  exit 0
fi

append_status "${ogname}" "OK" "${n_ids}" "${n_found}" "${n_after}" "ok"
echo "OK ${ogname} -> ${phy_out}"
