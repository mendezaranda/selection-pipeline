#!/bin/bash
#
#SBATCH --job-name=evigene_trinity
#SBATCH -N1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --gres=localtmp:80G
#SBATCH -t 12:00:00
#SBATCH -o /home/dmendez/outdir/%j_%x
#SBATCH -e /home/dmendez/errdir/%j_%x
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daniel.mendezaranda@mdc-berlin.de

# =============================================================================
# De novo transcriptome cleanup pipeline for Orthofinder + CDS downstream
# Daniel Mendez Aranda, 2026-02
#
# Input merged transcript FASTA may contain:
#   Curated:  >species:PR;...;NP_...;productName:GENE
#   Trinity:  >c11_g1_i2   or  >A_c..._g..._i... etc (may be multiple assemblies)
#
# Outputs:
#   outdir/<species>.primary_by_gene.pep.fa        (UNPREFIXED IDs)
#   outdir/<species>.primary_by_gene.cds.fa        (UNPREFIXED IDs; matches pep IDs)
#   outdir/<ORTHO_PREFIX>.orthofinder.pep.fa       (PREFIXED IDs)
#   outdir/<ORTHO_PREFIX>.orthofinder.cds.fa       (PREFIXED IDs; matches prefixed pep)
#   ORTHOFINDER_DIR/<ORTHO_PREFIX>.pep.fasta       (PREFIXED pep for OrthoFinder)
#   outdir/<species>.choice.tsv                    (trace)
#
# Important:
# - Curated TransDecoder.Predict can be overly strict -> we use LongOrfs and pick
#   one ORF per curated transcript, and also extract matching CDS.
# - Trinity IDs don't survive Evigene -> we collapse by TransDecoder GENE.<gene_id>
#   stripping t1/t2 to define locus.
#
# Optional de-dup across curated+trinity without DIAMOND:
#   DO_PROT_DEDUP=1 (default) + DEDUP_PROT_ID=0.99
#   This clusters merged proteins and keeps curated reps preferentially because curated
#   sequences are written first. CDS is filtered to kept IDs.
#
# Useful env vars:
#   ORTHO_PREFIX, BASE_OUTDIR, ORTHOFINDER_DIR, EVIGENEHOME
#   CDHIT_ID (cd-hit-est nucleotide), MIN_AA, COMPLETE_ONLY, KEEP_NONAME
#   DO_DIAMOND, DROP_TRINITY_MATCHED
#   DO_PROT_DEDUP=1, DEDUP_PROT_ID=0.99
# =============================================================================

set -eo pipefail

species="${1:-}"
input_fasta="${2:-}"
if [[ -z "$species" || -z "$input_fasta" ]]; then
  echo "Usage: $0 <species> <merged_transcripts_fasta>"
  exit 1
fi
if [[ ! -f "$input_fasta" ]]; then
  echo "ERROR: input fasta not found: $input_fasta"
  exit 1
fi

# Defaults (override via --export)
BASE_OUTDIR="${BASE_OUTDIR:-/fast/AG_Lewin/dmendez/transcriptomes/AMR/20260210}"
ORTHOFINDER_DIR="${ORTHOFINDER_DIR:-/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/01_proteomes}"
export EVIGENEHOME="${EVIGENEHOME:-/fast/AG_Lewin/dmendez/tools/evidentialgene}"

# Params/resources
NCPU="${SLURM_CPUS_PER_TASK:-16}"
CDHIT_ID="${CDHIT_ID:-0.97}"

MIN_AA="${MIN_AA:-100}"
COMPLETE_ONLY="${COMPLETE_ONLY:-0}"
KEEP_NONAME="${KEEP_NONAME:-1}"

DO_DIAMOND="${DO_DIAMOND:-0}"
DMND_MIN_BITSCORE="${DMND_MIN_BITSCORE:-60}"
DMND_MAX_TARGETS="${DMND_MAX_TARGETS:-1}"
DROP_TRINITY_MATCHED="${DROP_TRINITY_MATCHED:-0}"

DO_PROT_DEDUP="${DO_PROT_DEDUP:-1}"
DEDUP_PROT_ID="${DEDUP_PROT_ID:-0.99}"

ORTHO_PREFIX="${ORTHO_PREFIX:-$species}"

# Environment
source /home/dmendez/.bashrc
conda activate /fast/AG_Lewin/dmendez/.conda/envs/Trinity
module load mpi >/dev/null 2>&1 || true

if [[ -z "${EVIGENEHOME:-}" || ! -x "$EVIGENEHOME/scripts/prot/tr2aacds.pl" ]]; then
  echo "ERROR: Evigene not found at EVIGENEHOME=$EVIGENEHOME (missing tr2aacds.pl)"
  exit 1
fi

outdir="${BASE_OUTDIR}/${species}"
mkdir -p "${outdir}" "${outdir}/logs" "${ORTHOFINDER_DIR}"
cd "${outdir}"

# Files
merged_fa="${outdir}/${species}_merged.fa"
curated_tx="${outdir}/${species}.curated.transcripts.fa"
trinity_tx="${outdir}/${species}.trinity.transcripts.fa"

# Curated longorfs-derived
curated_longpep="${curated_tx}.transdecoder_dir/longest_orfs.pep"
curated_longcds="${curated_tx}.transdecoder_dir/longest_orfs.cds"
curated_td_pep="${outdir}/${species}.curated.longest_orfs.1perTx.pep"
curated_td_cds="${outdir}/${species}.curated.longest_orfs.1perTx.cds"
curated_orfmap="${outdir}/${species}.curated.tx_to_orf.tsv"

# Trinity (Evigene + TransDecoder Predict)
cdhit_fa="${outdir}/${species}_trinity_cdhit.uniq.fa"
okay_mrna="${outdir}/${species}_trinity_okay.mrna"
trinity_td_pep="${outdir}/${species}.trinity.transdecoder.pep"
trinity_td_cds="${outdir}/${species}.trinity.transdecoder.cds"

# Optional DIAMOND
dmnd_db="${outdir}/${species}.curated.dmnd"
dmnd_tsv="${outdir}/${species}.trinity_vs_curated.tsv"

# Final merged (unprefixed) + trace
merged_pep="${outdir}/${species}.primary_by_gene.pep.fa"
merged_cds="${outdir}/${species}.primary_by_gene.cds.fa"
pick_map="${outdir}/${species}.choice.tsv"

# Optional de-dup outputs
dedup_pep="${outdir}/${species}.primary_by_gene.dedup${DEDUP_PROT_ID}.pep.fa"
dedup_cds="${outdir}/${species}.primary_by_gene.dedup${DEDUP_PROT_ID}.cds.fa"
dedup_clstr="${dedup_pep}.clstr"
keep_ids="${outdir}/${species}.dedup.keep_ids.txt"

# Prefixed for orthofinder + prefixed cds
pref_pep="${outdir}/${ORTHO_PREFIX}.orthofinder.pep.fa"
pref_cds="${outdir}/${ORTHO_PREFIX}.orthofinder.cds.fa"
orthofinder_pep="${ORTHOFINDER_DIR}/${ORTHO_PREFIX}.pep.fasta"

count_fa () { grep -c '^>' "$1" 2>/dev/null || true; }

echo ">>> START"
echo "species: $species"
echo "ORTHO_PREFIX: $ORTHO_PREFIX"
echo "input: $input_fasta"
echo "outdir: $outdir"
echo "CDHIT_ID: $CDHIT_ID"
echo "MIN_AA: $MIN_AA"
echo "COMPLETE_ONLY: $COMPLETE_ONLY"
echo "KEEP_NONAME: $KEEP_NONAME"
echo "DO_DIAMOND: $DO_DIAMOND"
echo "DO_PROT_DEDUP: $DO_PROT_DEDUP (DEDUP_PROT_ID=$DEDUP_PROT_ID)"
date
echo

# -----------------------------------------------------------------------------
# Step 1: link merged
# -----------------------------------------------------------------------------
if [[ ! -s "${merged_fa}" ]]; then
  ln -sf "${input_fasta}" "${merged_fa}"
fi

# -----------------------------------------------------------------------------
# Step 2: split merged into curated vs trinity BEFORE Evigene
# -----------------------------------------------------------------------------
python3 - <<PY
from pathlib import Path

inp = Path("${merged_fa}")
cur = Path("${curated_tx}")
tri = Path("${trinity_tx}")

def it_fa(p):
    h=None; s=[]
    with p.open() as fh:
        for line in fh:
            line=line.rstrip("\n")
            if line.startswith(">"):
                if h: yield h, "".join(s)
                h=line.strip(); s=[]
            else:
                s.append(line.strip())
        if h: yield h, "".join(s)

ncur=ntri=0
with cur.open("w") as oc, tri.open("w") as ot:
    for h,seq in it_fa(inp):
        tok=h[1:].split()[0]
        if tok.startswith("species:PR;"):
            oc.write(h.split()[0]+"\n")
            for i in range(0,len(seq),60): oc.write(seq[i:i+60]+"\n")
            ncur+=1
        else:
            ot.write(h.split()[0]+"\n")
            for i in range(0,len(seq),60): ot.write(seq[i:i+60]+"\n")
            ntri+=1

print("curated transcripts:", ncur, "->", cur)
print("trinity/other transcripts:", ntri, "->", tri)
PY

# -----------------------------------------------------------------------------
# Step 3: Curated ORFs via TransDecoder.LongOrfs + pick 1 ORF per transcript (pep+cds)
# -----------------------------------------------------------------------------
if [[ -s "${curated_tx}" ]]; then
  echo ">>> TransDecoder.LongOrfs (curated)"
  TransDecoder.LongOrfs -t "${curated_tx}" 2>&1 | tee "${outdir}/logs/curated_longorfs.log"

  if [[ ! -s "${curated_longpep}" || ! -s "${curated_longcds}" ]]; then
    echo "ERROR: expected LongOrfs outputs missing: ${curated_longpep} and/or ${curated_longcds}"
    exit 1
  fi

  echo ">>> Curated: select 1 ORF per transcript (pep+cds)"
  python3 - <<PY
import re
from pathlib import Path

pep_in=Path("${curated_longpep}")
cds_in=Path("${curated_longcds}")
pep_out=Path("${curated_td_pep}")
cds_out=Path("${curated_td_cds}")
map_out=Path("${curated_orfmap}")

def fasta_iter(p):
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

def tx_key(tok):
    # remove TD suffixes conservatively, keep original curated transcript token intact
    tok=re.sub(r'\.p\d+$','',tok)
    tok=re.sub(r'(\.orf\d+|_orf\d+|\.m\.\d+)$','',tok)
    return tok

# pick best ORF (longest pep) per transcript key; store chosen ORF token
best={}  # tx -> (pep_len, orf_tok, pep_seq)
for h,seq in fasta_iter(pep_in):
    orf_tok=h[1:].split()[0]
    pep=seq.replace("*","")
    tx=tx_key(orf_tok)
    L=len(pep)
    prev=best.get(tx)
    if prev is None or L>prev[0]:
        best[tx]=(L, orf_tok, pep)

# load cds for chosen ORFs
chosen_orfs=set(v[1] for v in best.values())
cds_by_orf={}
for h,seq in fasta_iter(cds_in):
    orf_tok=h[1:].split()[0]
    if orf_tok in chosen_orfs:
        cds_by_orf[orf_tok]=seq

missing=[v[1] for v in best.values() if v[1] not in cds_by_orf]
if missing:
    # fail hard: pep/cds mismatch should not happen
    raise SystemExit(f"Missing CDS for {len(missing)} chosen curated ORFs; example: {missing[0]}")

with pep_out.open("w") as op, cds_out.open("w") as oc, map_out.open("w") as om:
    om.write("curated_tx\tchosen_orf\n")
    for tx,(L,orf_tok,pep) in best.items():
        # output headers as curated transcript token (tx), keep ORF token in map
        op.write(f">{tx}\n")
        for i in range(0,len(pep),60): op.write(pep[i:i+60]+"\n")
        cds=cds_by_orf[orf_tok]
        oc.write(f">{tx}\n")
        for i in range(0,len(cds),60): oc.write(cds[i:i+60]+"\n")
        om.write(f"{tx}\t{orf_tok}\n")

print("Curated transcripts with ORFs:", len(best))
print("Wrote pep:", pep_out)
print("Wrote cds:", cds_out)
print("Wrote map:", map_out)
PY
fi

# -----------------------------------------------------------------------------
# Step 4: Trinity cd-hit-est
# -----------------------------------------------------------------------------
if [[ -s "${trinity_tx}" && ! -s "${cdhit_fa}" ]]; then
  echo ">>> cd-hit-est (trinity/other)"
  cd-hit-est -i "${trinity_tx}" -o "${cdhit_fa}" \
    -c "${CDHIT_ID}" -aS 0.90 -g 1 -T "${NCPU}" -M 0 2>&1 | tee "${outdir}/cdhit.log"
fi

# -----------------------------------------------------------------------------
# Step 5: Trinity Evigene tr2aacds
# -----------------------------------------------------------------------------
if [[ -s "${cdhit_fa}" && ! -s "${okay_mrna}" ]]; then
  echo ">>> Evigene tr2aacds (trinity/other)"
  base="${species}_trinity_cdhit"
  evigene_work="${outdir}/${base}_evigene"
  mkdir -p "${evigene_work}"
  cd "${evigene_work}"
  ln -sf "${cdhit_fa}" "${base}.mrna"
  "$EVIGENEHOME/scripts/prot/tr2aacds.pl" \
    -mrna "${base}.mrna" \
    -NCPU "${NCPU}" \
    -MINCDS=90 2>&1 | tee "${outdir}/tr2aacds.log"

  ok_mrna_path="okayset/${base}.okay.mrna"
  if [[ -s "${ok_mrna_path}" ]]; then
    cp "${ok_mrna_path}" "${okay_mrna}"
  else
    echo "WARNING: Evigene produced no okay.mrna; falling back to cd-hit transcripts."
    cp "${cdhit_fa}" "${okay_mrna}"
  fi
  cd "${outdir}"
fi

# -----------------------------------------------------------------------------
# Step 6: Trinity TransDecoder Predict (pep+cds)
# -----------------------------------------------------------------------------
if [[ -s "${okay_mrna}" && ( ! -s "${trinity_td_pep}" || ! -s "${trinity_td_cds}" ) ]]; then
  echo ">>> TransDecoder (trinity/other)"
  TransDecoder.LongOrfs -t "${okay_mrna}" 2>&1 | tee "${outdir}/logs/trinity_longorfs.log"
  TransDecoder.Predict  -t "${okay_mrna}" 2>&1 | tee "${outdir}/logs/trinity_predict.log"
  cp -f "${okay_mrna}.transdecoder.pep" "${trinity_td_pep}"
  cp -f "${okay_mrna}.transdecoder.cds" "${trinity_td_cds}"
fi

# -----------------------------------------------------------------------------
# Step 7 (optional): DIAMOND trinity vs curated
# -----------------------------------------------------------------------------
if [[ "${DO_DIAMOND}" -eq 1 && -s "${curated_td_pep}" && -s "${trinity_td_pep}" ]]; then
  echo ">>> DIAMOND blastp (trinity vs curated)"
  conda activate orthofinder 2>/dev/null || true
  if command -v diamond >/dev/null 2>&1; then
    diamond makedb --in "${curated_td_pep}" -d "${dmnd_db}"
    diamond blastp -q "${trinity_td_pep}" -d "${dmnd_db}" \
      -o "${dmnd_tsv}" -f 6 qseqid sseqid bitscore pident length evalue \
      --max-target-seqs "${DMND_MAX_TARGETS}" --threads "${NCPU}"
  else
    echo "DIAMOND not found; skipping DIAMOND."
    DO_DIAMOND=0
  fi
  conda activate /fast/AG_Lewin/dmendez/.conda/envs/Trinity 2>/dev/null || true
fi

# -----------------------------------------------------------------------------
# Step 8: Build final primary proteins AND CDS (IDs match)
#   Curated IDs: collision-proof using productName + NP/comp/Contig or hash
#   Trinity IDs: locus from TD header GENE.<gene_id> stripping t1/t2
# -----------------------------------------------------------------------------
python3 - <<PY
import re, hashlib
from pathlib import Path
from collections import defaultdict

cur_pep = Path("${curated_td_pep}")
cur_cds = Path("${curated_td_cds}")
tri_pep = Path("${trinity_td_pep}")
tri_cds = Path("${trinity_td_cds}")
dmnd_tsv = Path("${dmnd_tsv}")
pep_out = Path("${merged_pep}")
cds_out = Path("${merged_cds}")
pick = Path("${pick_map}")

min_aa = int("${MIN_AA}")
complete_only = int("${COMPLETE_ONLY}")
keep_noname = int("${KEEP_NONAME}")
do_dmnd = int("${DO_DIAMOND}")
min_bitscore = float("${DMND_MIN_BITSCORE}")
drop_trinity_matched = int("${DROP_TRINITY_MATCHED}")

hdr_gene = re.compile(r'^>(\S+).*?GENE\.([^~\s]+)~~([^\s]+)', re.IGNORECASE)
hdr_scorelen = re.compile(r'score=([0-9.]+).*?len:?[\s]*([0-9]+)', re.IGNORECASE)
orf_type_re = re.compile(r'ORF\s+type:([^\s]+)', re.IGNORECASE)

acc_re = re.compile(r'(NP_\d+(?:\.\d+)?)')
prod_re = re.compile(r'productName:([^;\s]+)')
comp_re = re.compile(r'(comp\d+_seq\d+|Contig\d+)', re.IGNORECASE)

def fasta_iter(path):
    if not path.exists() or path.stat().st_size==0:
        return
    h=None; s=[]
    with path.open() as fh:
        for line in fh:
            line=line.rstrip("\n")
            if line.startswith(">"):
                if h: yield h,"".join(s)
                h=line; s=[]
            else:
                s.append(line.strip())
        if h: yield h,"".join(s)

def load_fa_dict(path):
    d={}
    for h,seq in fasta_iter(path):
        d[h[1:].split()[0]] = seq
    return d

def clean(x:str)->str:
    x=x.strip()
    x=re.sub(r'[^A-Za-z0-9_.-]+','_',x)
    x=re.sub(r'_+','_',x).strip('_')
    return x or "UNKNOWN"

def meta(header, seq):
    seq = seq.replace("*","")
    L = len(seq)
    score = 0.0
    ms = hdr_scorelen.search(header)
    if ms:
        try: score=float(ms.group(1))
        except: score=0.0
        try: L=int(ms.group(2))
        except: L=len(seq)
    ot = orf_type_re.search(header)
    if ot:
        return score, L, ot.group(1).lower(), True
    return score, L, "unknown", False

def curated_uid(cur_token: str) -> str:
    pn = prod_re.search(cur_token)
    np = acc_re.search(cur_token)
    cp = comp_re.search(cur_token)
    if pn and np:
        return clean(f"{pn.group(1)}__{np.group(1)}")
    if pn and cp:
        return clean(f"{pn.group(1)}__{cp.group(1)}")
    if pn:
        h = hashlib.md5(cur_token.encode()).hexdigest()[:10]
        return clean(f"{pn.group(1)}__{h}")
    if np:
        return clean(np.group(1))
    h = hashlib.md5(cur_token.encode()).hexdigest()[:10]
    return clean(f"cur__{h}")

# Load CDS dictionaries
cur_cds_d = load_fa_dict(cur_cds) if cur_cds.exists() else {}
tri_cds_d = load_fa_dict(tri_cds) if tri_cds.exists() else {}

# Curated: pick best per curated_uid (longest pep; tie by score if possible)
cur_best = {}  # uid -> (L, score, src_token, pep_seq, cds_seq)
cur_stats = {"total":0,"minlen_skip":0,"complete_skip":0,"no_cds":0}

for h,pep in fasta_iter(cur_pep):
    cur_stats["total"] += 1
    tok=h[1:].split()[0]   # curated transcript token
    pep=pep.replace("*","")
    if min_aa and len(pep)<min_aa:
        cur_stats["minlen_skip"] += 1
        continue

    score,L,orf_type,has_orf = meta(h, pep)
    if complete_only and has_orf and orf_type != "complete":
        cur_stats["complete_skip"] += 1
        continue

    cds = cur_cds_d.get(tok)
    if cds is None:
        cur_stats["no_cds"] += 1
        continue

    uid = curated_uid(tok)
    if (not keep_noname) and re.match(r'^NonamEVm\d+', uid):
        continue

    prev = cur_best.get(uid)
    if prev is None or (L, score) > (prev[0], prev[1]):
        cur_best[uid] = (L, score, tok, pep, cds)

# Trinity: group by locus from TD gene_id, choose representative
tri_by_locus = defaultdict(list)  # locus -> list[(full_id, gene_id, L, score, pep, cds)]
tri_stats = {"total":0,"minlen_skip":0,"complete_skip":0,"no_geneid":0,"no_cds":0}

for h,pep in fasta_iter(tri_pep):
    tri_stats["total"] += 1
    full_tok=h[1:].split()[0]
    pep=pep.replace("*","")
    if min_aa and len(pep)<min_aa:
        tri_stats["minlen_skip"] += 1
        continue

    score,L,orf_type,has_orf = meta(h, pep)
    if complete_only and has_orf and orf_type != "complete":
        tri_stats["complete_skip"] += 1
        continue

    mg = hdr_gene.match(h)
    if mg:
        gene_id = mg.group(2)
        locus = re.sub(r't\d+$','', gene_id)
    else:
        tri_stats["no_geneid"] += 1
        gene_id = "NA"
        locus = full_tok

    if (not keep_noname) and re.match(r'^NonamEVm\d+', locus):
        continue

    cds = tri_cds_d.get(full_tok)
    if cds is None:
        tri_stats["no_cds"] += 1
        continue

    tri_by_locus[locus].append((full_tok, gene_id, L, score, pep, cds))

# Optional DIAMOND: choose best q per locus
best_q={}
best_bs={}
if do_dmnd and dmnd_tsv.exists() and dmnd_tsv.stat().st_size>0:
    q2locus={}
    for locus,lst in tri_by_locus.items():
        for full_tok,gene_id,L,score,pep,cds in lst:
            q2locus[full_tok]=locus
    with dmnd_tsv.open() as fh:
        for line in fh:
            q,s,bs,pi,aln,ev = line.rstrip("\n").split("\t")
            bs=float(bs)
            if bs < min_bitscore:
                continue
            locus=q2locus.get(q)
            if locus is None:
                continue
            if locus not in best_bs or bs>best_bs[locus]:
                best_bs[locus]=bs
                best_q[locus]=q

tri_chosen={}  # locus -> (full_tok,gene_id,L,score,reason,bitscore,pep,cds)
for locus,lst in tri_by_locus.items():
    if locus in best_q:
        q=best_q[locus]
        hit=None
        for full_tok,gene_id,L,score,pep,cds in lst:
            if full_tok==q:
                hit=(full_tok,gene_id,L,score,"match_curated",best_bs.get(locus,0.0),pep,cds)
                break
        if hit is None:
            full_tok,gene_id,L,score,pep,cds=max(lst, key=lambda x:(x[2],x[3]))
            hit=(full_tok,gene_id,L,score,"longest_fallback",0.0,pep,cds)
        tri_chosen[locus]=hit
    else:
        full_tok,gene_id,L,score,pep,cds=max(lst, key=lambda x:(x[2],x[3]))
        tri_chosen[locus]=(full_tok,gene_id,L,score,"longest",0.0,pep,cds)

# Write merged pep+cds with identical IDs
used=set()
def unique(n):
    if n not in used:
        used.add(n); return n
    k=2
    while f"{n}_dup{k}" in used:
        k+=1
    nn=f"{n}_dup{k}"
    used.add(nn)
    return nn

pep_out.parent.mkdir(parents=True, exist_ok=True)
with pep_out.open("w") as op, cds_out.open("w") as oc:
    # curated first
    for uid,(L,score,src,pep,cds) in sorted(cur_best.items(), key=lambda x:(-x[1][0],-x[1][1],x[0])):
        name=unique(uid)
        op.write(f">{name}\n")
        for i in range(0,len(pep),60): op.write(pep[i:i+60]+"\n")
        oc.write(f">{name}\n")
        for i in range(0,len(cds),60): oc.write(cds[i:i+60]+"\n")

    # trinity
    for locus,(full_tok,gene_id,L,score,reason,bitscore,pep,cds) in sorted(tri_chosen.items(), key=lambda x:(-x[1][2],-x[1][3],x[0])):
        if drop_trinity_matched and reason=="match_curated":
            continue
        name=unique(clean(locus))
        op.write(f">{name}\n")
        for i in range(0,len(pep),60): op.write(pep[i:i+60]+"\n")
        oc.write(f">{name}\n")
        for i in range(0,len(cds),60): oc.write(cds[i:i+60]+"\n")

with pick.open("w") as ph:
    ph.write("type\tfinal_id\tsrc_full_id\tgene_id\tlen\tscore\treason\tbitscore\n")
    for uid,(L,score,src,pep,cds) in sorted(cur_best.items()):
        ph.write(f"curated\t{uid}\t{src}\tNA\t{L}\t{score}\tkept\t0\n")
    for locus,(full_tok,gene_id,L,score,reason,bitscore,pep,cds) in sorted(tri_chosen.items()):
        ph.write(f"trinity\t{clean(locus)}\t{full_tok}\t{gene_id}\t{L}\t{score}\t{reason}\t{bitscore}\n")

print("Curated stats:", cur_stats, "unique_uids:", len(cur_best))
print("Trinity stats:", tri_stats, "unique_loci:", len(tri_chosen))
print("Wrote pep:", pep_out)
print("Wrote cds:", cds_out)
print("Trace:", pick)
PY

# -----------------------------------------------------------------------------
# Step 9: Optional protein-level de-dup (cd-hit) + filter CDS to kept IDs
# -----------------------------------------------------------------------------
final_pep="${merged_pep}"
final_cds="${merged_cds}"

if [[ "${DO_PROT_DEDUP}" -eq 1 ]]; then
  echo ">>> Protein de-dup with cd-hit (c=${DEDUP_PROT_ID})"
  if ! command -v cd-hit >/dev/null 2>&1; then
    echo "ERROR: cd-hit not found in PATH (need cd-hit for protein de-dup)."
    exit 1
  fi

  cd-hit -i "${merged_pep}" -o "${dedup_pep}" -c "${DEDUP_PROT_ID}" -n 5 -T "${NCPU}" -M 0 \
    2>&1 | tee "${outdir}/logs/cdhit_protein_dedup.log"

  if [[ ! -s "${dedup_clstr}" ]]; then
    echo "ERROR: missing cd-hit cluster file: ${dedup_clstr}"
    exit 1
  fi

  echo ">>> Extract kept IDs from cd-hit clusters"
  python3 - <<PY
from pathlib import Path
cl=Path("${dedup_clstr}")
out=Path("${keep_ids}")

kept=[]
current=[]
for line in cl.read_text().splitlines():
    if line.startswith(">Cluster"):
        if current:
            rep=[x for x in current if x[1]]
            if rep:
                kept.append(rep[0][0])
            else:
                kept.append(current[0][0])
        current=[]
        continue
    # example: 0 123aa, >ID... *
    line=line.strip()
    if not line:
        continue
    m=line.split(">")
    if len(m)<2:
        continue
    rest=m[1]
    pid=rest.split("...")[0].strip()
    is_rep=line.endswith("*")
    current.append((pid,is_rep))

# last cluster
if current:
    rep=[x for x in current if x[1]]
    kept.append(rep[0][0] if rep else current[0][0])

out.write_text("\n".join(kept) + "\n")
print("Kept IDs:", len(kept), "->", out)
PY

  echo ">>> Filter CDS to kept IDs"
  python3 - <<PY
from pathlib import Path

ids=set(Path("${keep_ids}").read_text().split())
inp=Path("${merged_cds}")
out=Path("${dedup_cds}")

def fasta_iter(p):
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

with out.open("w") as oh:
    for h,seq in fasta_iter(inp):
        pid=h[1:].split()[0]
        if pid in ids:
            oh.write(f">{pid}\n")
            for i in range(0,len(seq),60):
                oh.write(seq[i:i+60]+"\n")

print("Wrote filtered CDS:", out)
PY

  final_pep="${dedup_pep}"
  final_cds="${dedup_cds}"
fi

# -----------------------------------------------------------------------------
# Step 10: Prefix IDs for OrthoFinder and copy pep; also write prefixed CDS
# -----------------------------------------------------------------------------
awk -v sp="${ORTHO_PREFIX}" '
  /^>/ { id=substr($0,2); sub(/[[:space:]].*$/, "", id); print ">" sp "_" id; next }
  { print }
' "${final_pep}" > "${pref_pep}"

awk -v sp="${ORTHO_PREFIX}" '
  /^>/ { id=substr($0,2); sub(/[[:space:]].*$/, "", id); print ">" sp "_" id; next }
  { print }
' "${final_cds}" > "${pref_cds}"

cp "${pref_pep}" "${orthofinder_pep}"

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------
echo
echo "Counts:"
echo "merged transcripts     : $(count_fa "${merged_fa}")"
echo "curated transcripts    : $(count_fa "${curated_tx}")"
echo "trinity transcripts    : $(count_fa "${trinity_tx}")"
echo "curated 1perTx pep     : $(count_fa "${curated_td_pep}")"
echo "curated 1perTx cds     : $(count_fa "${curated_td_cds}")"
echo "trinity TD pep         : $(count_fa "${trinity_td_pep}")"
echo "trinity TD cds         : $(count_fa "${trinity_td_cds}")"
echo "final primary pep      : $(count_fa "${final_pep}")"
echo "final primary cds      : $(count_fa "${final_cds}")"
echo "prefixed pep           : $(count_fa "${pref_pep}")"
echo "prefixed cds           : $(count_fa "${pref_cds}")"
echo
echo "Primary (pep):   ${final_pep}"
echo "Primary (cds):   ${final_cds}"
echo "OrthoFinder pep: ${orthofinder_pep}"
echo "OrthoFinder cds: ${pref_cds}"
echo "Trace:          ${pick_map}"
echo "Curated ORF map:${curated_orfmap}"
echo
echo ">>> DONE"
date
