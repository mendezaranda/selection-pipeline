#!/usr/bin/env bash
set -eo pipefail
# If configs/env.sh has already been sourced in your shell (see
# configs/env.example.sh), CODEML/TRIMAL below pick that up automatically.
# Otherwise export them yourself, or rely on the PATH-based defaults.

# Usage:
#   bash og_to_codon_paml_codeml.sh <OGxxxx.fa> <ALL_CDS.fa> <OUTDIR> <TREE_FG.nwk> [THREADS]
OG_FA="${1:?Usage: bash og_to_codon_paml_codeml.sh <OG.fa> <ALL_CDS.fa> <OUTDIR> <TREE_FG.nwk> [THREADS]}"
ALL_CDS="${2:?Usage: bash og_to_codon_paml_codeml.sh <OG.fa> <ALL_CDS.fa> <OUTDIR> <TREE_FG.nwk> [THREADS]}"
OUTDIR="${3:?Usage: bash og_to_codon_paml_codeml.sh <OG.fa> <ALL_CDS.fa> <OUTDIR> <TREE_FG.nwk> [THREADS]}"
TREE_FG="${4:?Usage: bash og_to_codon_paml_codeml.sh <OG.fa> <ALL_CDS.fa> <OUTDIR> <TREE_FG.nwk> [THREADS]}"
THREADS="${5:-8}"

# Tunables (override via env, or set in configs/env.sh)
MIN_TAXA="${MIN_TAXA:-10}"
DO_TRIMAL="${DO_TRIMAL:-1}"
TRIMAL="${TRIMAL:-trimal}"
FILTER_STOPS="${FILTER_STOPS:-1}"
FILTER_FRAME="${FILTER_FRAME:-1}"

# Run codeml?
RUN_CODEML="${RUN_CODEML:-1}"     # 1 run codeml, 0 just write alignments
CODEML="${CODEML:-codeml}"

command -v mafft >/dev/null 2>&1 || { echo "mafft not found"; exit 2; }
if [ "${DO_TRIMAL}" -eq 1 ]; then
  command -v "${TRIMAL}" >/dev/null 2>&1 || { echo "trimAl not found at ${TRIMAL}"; exit 2; }
fi
if [ "${RUN_CODEML}" -eq 1 ]; then
  command -v "${CODEML}" >/dev/null 2>&1 || { echo "codeml not found in PATH (CODEML=${CODEML})"; exit 2; }
fi
[ -s "${TREE_FG}" ] || { echo "Missing tree file: ${TREE_FG}"; exit 2; }

mkdir -p "${OUTDIR}/tmp" "${OUTDIR}/aa_aln" "${OUTDIR}/codon_aln" "${OUTDIR}/paml_phy" "${OUTDIR}/hyphy_fa" "${OUTDIR}/logs" "${OUTDIR}/codeml_runs"

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
trimal_cols="${work}/${ogname}.trimal_cols.txt"
codon_fa="${OUTDIR}/codon_aln/${ogname}.codon.aln.fa"

paml_phy="${OUTDIR}/paml_phy/${ogname}.codon.phy"
hyphy_fa="${OUTDIR}/hyphy_fa/${ogname}.codon.fasta"
log="${OUTDIR}/logs/${ogname}.log"

# Per-task codeml summary (no locking)
TASK_TAG="${SLURM_ARRAY_TASK_ID:-0}"
SUMMARY_TSV="${OUTDIR}/codeml_summary.task${TASK_TAG}.tsv"

append_summary () {
  local og="$1" st="$2" ntaxa="$3" nsites="$4" ln0="$5" ln1="$6" lrt="$7" p="$8" note="$9"
  if [ ! -s "${SUMMARY_TSV}" ]; then
    echo -e "OG\tstatus\tntaxa\tnsites\tlnL_null\tlnL_alt\tLRT\tp_df1\tnote" > "${SUMMARY_TSV}"
  fi
  echo -e "${og}\t${st}\t${ntaxa}\t${nsites}\t${ln0}\t${ln1}\t${lrt}\t${p}\t${note}" >> "${SUMMARY_TSV}"
}

# PAML's lnL line looks like:
#   lnL(ntime: 43  np: 45):  -12345.678901      +0.000000
# i.e. the value follows "):" -- there is no "=" on this line in standard
# PAML 4.x output. Used both to decide whether a previous codeml run for this
# OG actually succeeded (skip-check below) and to parse the final result.
parse_lnl() {
  local f="$1" v
  v="$(grep -m1 -oE 'lnL\([^)]*\):[[:space:]]*-?[0-9]+\.[0-9]+' "$f" 2>/dev/null \
       | grep -oE -- '-?[0-9]+\.[0-9]+$')"
  if [ -z "${v}" ]; then
    # Fallback for older/alternate builds that do print "lnL = <value>"
    v="$(awk '/lnL/ && /=/{for(i=1;i<=NF;i++){if($i=="="){print $(i+1); exit}}}' "$f" 2>/dev/null | head -n1)"
  fi
  echo "${v}"
}

# Skip if already done. NOTE: a non-empty codeml_alt.out is NOT enough to call
# an OG done -- codeml writes codon-usage tables and distance matrices before
# it ever reaches the ML tree/lnL step, so a run that aborted partway through
# (e.g. the tree/taxon-set mismatch bug 3 fixed below) still leaves a
# non-empty, but incomplete, codeml_alt.out. Require an actual parseable lnL
# so a rerun after a bugfix reprocesses exactly the OGs that never really
# finished, and leaves already-correct OGs alone (cheap incremental reruns).
codeml_alt_out="${OUTDIR}/codeml_runs/${ogname}/codeml_alt.out"
codeml_ok=0
if [ -s "${codeml_alt_out}" ] && [ -n "$(parse_lnl "${codeml_alt_out}")" ]; then
  codeml_ok=1
fi

if [ -s "${paml_phy}" ] && [ -s "${hyphy_fa}" ] && { [ "${RUN_CODEML}" -eq 0 ] || [ "${codeml_ok}" -eq 1 ] || [ -e "${OUTDIR}/codeml_runs/${ogname}/SKIPPED_NO_FOREGROUND" ]; }; then
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
  append_summary "${ogname}" "SKIP_TOO_FEW_CDS" "${n_cds}" "NA" "NA" "NA" "NA" "NA" "too_few_cds"
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
  append_summary "${ogname}" "SKIP_TOO_FEW_AFTER_FILTER" "${n_filt}" "NA" "NA" "NA" "NA" "NA" "too_few_after_filter"
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
    if aas and aas[-1] == "*":
        aas = aas[:-1]
    return "".join(aas)

with open(out,"w") as o:
    for name, seq in read_fasta(inp):
        seq=seq.upper().replace("U","T")
        o.write(f">{name}\n{translate(seq)}\n")
PY

# 5) Align AA, optional trim.
#    IMPORTANT: we keep BOTH the untrimmed alignment (aa_aln, which is in exact
#    1:1 codon correspondence with cds_filt) and, when trimming, the trimAl
#    "-colnumbering" map of which alignment columns survived. Backtranslation
#    (step 6) uses that column map instead of guessing from sequence length,
#    because trimAl's -automated1 can drop columns from the MIDDLE of the
#    alignment, not just the ends.
mafft --auto --thread "${THREADS}" "${aa_from_cds}" > "${aa_aln}"
: > "${trimal_cols}"
if [ "${DO_TRIMAL}" -eq 1 ]; then
  "${TRIMAL}" -in "${aa_aln}" -out "${aa_trim}" -automated1 -colnumbering > "${trimal_cols}" 2>&1
else
  cp "${aa_aln}" "${aa_trim}"
fi

# 6) Backtranslate to codon alignment and rename to species-only.
#    If an OG contains >1 sequence for the same species, skip that OG.
python3 - <<PY >> "${log}" 2>&1
import sys, re, os

AA_ALN_FULL="${aa_aln}"      # untrimmed; 1:1 codon correspondence with cds_filt
AA_ALN_TRIM="${aa_trim}"     # trimAl output (used only to sanity-check width)
CDS="${cds_filt}"
COLS_FILE="${trimal_cols}"
DO_TRIMAL=int("${DO_TRIMAL}")
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

full_aln = {n:s for n,s in read_fasta(AA_ALN_FULL)}
trim_aln = {n:s for n,s in read_fasta(AA_ALN_TRIM)}
cds = {n:s.upper().replace("U","T") for n,s in read_fasta(CDS)}

aln_width = len(next(iter(full_aln.values()))) if full_aln else 0

if DO_TRIMAL:
    text = open(COLS_FILE).read() if COLS_FILE and os.path.exists(COLS_FILE) else ""
    kept_cols = [int(x) for x in re.findall(r"\d+", text)]
else:
    kept_cols = list(range(aln_width))

trim_width = len(next(iter(trim_aln.values()))) if trim_aln else None
if trim_width is None or len(kept_cols) != trim_width or (kept_cols and max(kept_cols) >= aln_width):
    # Column map didn't parse cleanly against the trimmed alignment we actually
    # got back -- bail out for this OG rather than risk misaligned codons.
    open(OUT,"w").close()
    sys.stderr.write(
        f"SKIP {os.path.basename(OUT)}: trimAl column map unusable "
        f"(parsed {len(kept_cols)} cols, trimmed width {trim_width}, aln width {aln_width})\n"
    )
    sys.exit(0)

out={}
bad=[]
dup_species=set()

for n, aln in full_aln.items():
    if n not in cds:
        bad.append((n,"no_cds"))
        continue
    cds_seq = cds[n]

    # Walk the FULL (untrimmed) alignment column-by-column, pulling codons off
    # cds_seq in order for every non-gap position. This is exact because
    # aa_from_cds (and hence aa_aln) was translated directly from cds_filt.
    codon_of_col = []
    ci = 0
    ok = True
    for ch in aln:
        if ch in "-.":
            codon_of_col.append("---")
            continue
        if (ci+1)*3 > len(cds_seq):
            bad.append((n,"cds_too_short"))
            ok = False
            break
        c = cds_seq[ci*3:(ci+1)*3]
        a = aa_from_codon(c)
        if ch != "X" and a != "X" and ch != a:
            bad.append((n,f"aa_mismatch_col{len(codon_of_col)}"))
            ok = False
            break
        codon_of_col.append(c)
        ci += 1
    if not ok:
        continue

    try:
        trimmed_codons = [codon_of_col[c] for c in kept_cols]
    except IndexError:
        bad.append((n,"col_index_out_of_range"))
        continue

    seq_out = "".join(trimmed_codons)
    aa_check = "".join(aa_from_codon(seq_out[i:i+3]) for i in range(0,len(seq_out),3))
    if "*" in aa_check[:-1]:
        bad.append((n,"internal_stop_after_trim"))
        continue

    sp = n.split("_",1)[0]
    if sp in out:
        dup_species.add(sp)
        continue
    out[sp] = seq_out

if dup_species:
    open(OUT,"w").close()
    sys.stderr.write("Duplicate species in OG, skipping: " + ",".join(sorted(dup_species)) + "\n")
    sys.exit(0)

if len(out) < MIN_TAXA:
    open(OUT,"w").close()
    sys.exit(0)

L=None
for sp, seq in out.items():
    if L is None: L=len(seq)
    if len(seq)!=L:
        open(OUT,"w").close()
        sys.exit(0)

if L is not None and L % 3 != 0:
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
  append_summary "${ogname}" "SKIP_NO_CODON_ALN" "NA" "NA" "NA" "NA" "NA" "NA" "no_codon_aln"
  exit 0
fi

# 7) Write HyPhy FASTA + PAML PHYLIP
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

# 8) Run codeml per OG (one-ratio null vs two-ratio alt; foreground branches marked #1 in TREE_FG)
if [ "${RUN_CODEML}" -eq 1 ]; then
  ogdir="${OUTDIR}/codeml_runs/${ogname}"
  mkdir -p "${ogdir}"
  cp -f "${paml_phy}" "${ogdir}/alignment.phy"

  read -r ntaxa nsites < <(awk 'NR==1{print $1"\t"$2; exit}' "${ogdir}/alignment.phy")

  # Prune the full-species foreground tree (TREE_FG) down to exactly the taxa
  # present in THIS OG's alignment, keeping the PAML "#1" foreground mark on
  # whichever leaf/clade carries it. codeml requires an exact taxon-set match
  # between the tree and the alignment; not every OG has all species (missing
  # single-copy orthologs is normal), so reusing the same full tree for every
  # OG made codeml abort ("Error: check #seqs and tree: perhaps too many '('?")
  # for any OG missing even one species relative to TREE_FG.
  if python3 - "${codon_fa}" "${TREE_FG}" "${ogdir}/tree_fg.nwk" <<'PY' > "${ogdir}/prune_tree.log" 2>&1
import sys, re

CODON_FA, TREE_IN, TREE_OUT = sys.argv[1], sys.argv[2], sys.argv[3]

def read_fasta_names(fn):
    names = []
    for line in open(fn):
        if line.startswith(">"):
            names.append(line[1:].strip().split()[0])
    return names

keep = set(read_fasta_names(CODON_FA))

def parse_newick(s):
    s = s.strip()
    if s.endswith(';'):
        s = s[:-1]
    pos = [0]

    def parse_label():
        start = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ',();:':
            pos[0] += 1
        return s[start:pos[0]]

    def split_tag(label):
        m = re.match(r'^(.*?)(#\d+)?$', label)
        return m.group(1), m.group(2)

    def parse_node():
        node = {'name': None, 'tag': None, 'length': None, 'children': []}
        if s[pos[0]] == '(':
            pos[0] += 1
            node['children'].append(parse_node())
            while s[pos[0]] == ',':
                pos[0] += 1
                node['children'].append(parse_node())
            assert s[pos[0]] == ')', f"expected ) at {pos[0]}"
            pos[0] += 1
            label = parse_label()
            name, tag = split_tag(label)
            node['name'] = name if name else None
            node['tag'] = tag
        else:
            label = parse_label()
            name, tag = split_tag(label)
            node['name'] = name
            node['tag'] = tag
        if pos[0] < len(s) and s[pos[0]] == ':':
            pos[0] += 1
            node['length'] = parse_label()
        return node

    return parse_node()

def leaves(node, out=None):
    if out is None: out = []
    if not node['children']:
        out.append(node)
    else:
        for c in node['children']:
            leaves(c, out)
    return out

def prune(node, keep):
    if not node['children']:
        return dict(node) if node['name'] in keep else None
    new_children = []
    for c in node['children']:
        pc = prune(c, keep)
        if pc is not None:
            new_children.append(pc)
    if not new_children:
        return None
    if len(new_children) == 1:
        child = new_children[0]
        def add_len(a, b):
            fa = float(a) if a not in (None, '') else 0.0
            fb = float(b) if b not in (None, '') else 0.0
            return repr(fa + fb)
        merged_len = add_len(node.get('length'), child.get('length'))
        child = dict(child)
        child['length'] = merged_len
        return child
    return {'name': node['name'], 'tag': node['tag'], 'length': node['length'], 'children': new_children}

def to_newick(node):
    def tag_str(n):
        return n['tag'] if n['tag'] else ''
    def rec(n):
        if not n['children']:
            s = (n['name'] or '') + tag_str(n)
        else:
            s = '(' + ','.join(rec(c) for c in n['children']) + ')' + (n['name'] or '') + tag_str(n)
        if n.get('length') not in (None, ''):
            s += ':' + n['length']
        return s
    return rec(node) + ';'

tree = parse_newick(open(TREE_IN).read())
pruned = prune(tree, keep)
if pruned is None:
    print("NO_TAXA_LEFT")
    sys.exit(3)

pruned_leaves = leaves(pruned)
pruned_names = set(l['name'] for l in pruned_leaves)
if pruned_names != keep:
    print(f"TAXON_SET_MISMATCH: pruned={sorted(pruned_names)} keep={sorted(keep)}")
    sys.exit(4)

fg = [l for l in pruned_leaves if l['tag']]
if not fg:
    print("NO_FOREGROUND_TAXON")
    sys.exit(2)

with open(TREE_OUT, "w") as o:
    o.write(to_newick(pruned) + "\n")
print(f"OK: pruned to {len(pruned_names)} taxa, foreground={[l['name'] for l in fg]}")
PY
  then
    prune_rc=0
  else
    prune_rc=$?
  fi

  if [ "${prune_rc}" -eq 2 ]; then
    echo "SKIP ${ogname} (foreground taxon absent from this OG's alignment; see ${ogdir}/prune_tree.log)"
    touch "${ogdir}/SKIPPED_NO_FOREGROUND"
    append_summary "${ogname}" "SKIP_NO_FOREGROUND_TAXON" "${ntaxa}" "${nsites}" "NA" "NA" "NA" "NA" "foreground_taxon_missing"
  elif [ "${prune_rc}" -ne 0 ]; then
    echo "SKIP ${ogname} (tree pruning failed rc=${prune_rc}; see ${ogdir}/prune_tree.log)"
    append_summary "${ogname}" "SKIP_TREE_PRUNE_FAILED" "${ntaxa}" "${nsites}" "NA" "NA" "NA" "NA" "tree_prune_failed"
  else

  cat > "${ogdir}/codeml_null.ctl" <<'CTL'
      seqfile = alignment.phy
     treefile = tree_fg.nwk
      outfile = codeml_null.out

        noisy = 3
      verbose = 0
      runmode = 0

      seqtype = 1
    CodonFreq = 2

        clock = 0
       aaDist = 0

      model = 0
    NSsites = 0

  fix_kappa = 0
      kappa = 2

  fix_omega = 0
      omega = 0.2

    cleandata = 0
CTL

  cat > "${ogdir}/codeml_alt.ctl" <<'CTL'
      seqfile = alignment.phy
     treefile = tree_fg.nwk
      outfile = codeml_alt.out

        noisy = 3
      verbose = 0
      runmode = 0

      seqtype = 1
    CodonFreq = 2

        clock = 0
       aaDist = 0

      model = 2
    NSsites = 0

  fix_kappa = 0
      kappa = 2

  fix_omega = 0
      omega = 0.2

    cleandata = 0
CTL

  ( cd "${ogdir}" && "${CODEML}" codeml_null.ctl > codeml_null.log 2>&1 ) || true
  ( cd "${ogdir}" && "${CODEML}" codeml_alt.ctl  > codeml_alt.log  2>&1 ) || true

  # parse_lnl() is defined near the top of this script (also used by the
  # skip-check above).
  lnl_null="$(parse_lnl "${ogdir}/codeml_null.out" || true)"
  lnl_alt="$(parse_lnl "${ogdir}/codeml_alt.out"  || true)"

  if [ -n "${lnl_null}" ] && [ -n "${lnl_alt}" ]; then
    lrt="$(python3 - <<PY
ln0=float("${lnl_null}")
ln1=float("${lnl_alt}")
print(max(0.0, 2.0*(ln1-ln0)))
PY
)"
    p="$(python3 - <<PY
import math
x=float("${lrt}")
p = math.erfc(math.sqrt(x/2.0)) if x >= 0 else 1.0
print(p)
PY
)"
    append_summary "${ogname}" "OK" "${ntaxa}" "${nsites}" "${lnl_null}" "${lnl_alt}" "${lrt}" "${p}" "codeml_branch_model"
  else
    append_summary "${ogname}" "CODEML_FAIL" "${ntaxa}" "${nsites}" "${lnl_null:-NA}" "${lnl_alt:-NA}" "NA" "NA" "parse_lnl_failed"
  fi
  fi
fi

echo "OK ${ogname} -> ${paml_phy} ${hyphy_fa}"
