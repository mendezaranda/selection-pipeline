#!/bin/bash

set -eo pipefail

# ==============================================================================
# extract_mouse_gene_ogs_and_scan.sh
#
# Inputs:
#   1) ORTHO_RESULTS_DIR : OrthoFinder Results_* directory containing Orthogroups/ + sequences
#   2) MAPPING_TSV       : MusMus_protein2gene.tsv (protein_id<TAB>gene)
#   3) QUERY_LIST        : file with one entry per line: either NP_/XP_ or gene symbol
#   4) OUTDIR            : output directory
#
# Env (optional):
#   MOUSE_COL            : column name in Orthogroups.tsv for mouse (default: MusMus)
#   CLADE3               : comma-separated target trio (default: HelKap,HetGla,GeoCap)
#   FOCAL                : focal species for "private" changes (default: HelKap)
#
# Requirements:
#   - mafft on PATH (conda install -c bioconda mafft)
#
# Outputs (in OUTDIR):
#   - gene_orthogroups.tsv
#   - alignments/<GENE>__<OG>.aa.aln.fa
#   - reports/<GENE>__<OG>.mutations.tsv
#   - reports/summary.tsv
# ==============================================================================

ORTHO_RESULTS_DIR="${1:-}"
MAPPING_TSV="${2:-}"
QUERY_LIST="${3:-}"
OUTDIR="${4:-}"

if [[ -z "${ORTHO_RESULTS_DIR}" || -z "${MAPPING_TSV}" || -z "${QUERY_LIST}" || -z "${OUTDIR}" ]]; then
  echo "Usage: $0 <ORTHO_RESULTS_DIR> <MusMus_protein2gene.tsv> <query_list.txt> <outdir>"
  exit 1
fi

if [[ ! -d "${ORTHO_RESULTS_DIR}" ]]; then
  echo "ERROR: ORTHO_RESULTS_DIR not found: ${ORTHO_RESULTS_DIR}"
  exit 2
fi
if [[ ! -s "${MAPPING_TSV}" ]]; then
  echo "ERROR: mapping TSV not found/empty: ${MAPPING_TSV}"
  exit 2
fi
if [[ ! -s "${QUERY_LIST}" ]]; then
  echo "ERROR: query list not found/empty: ${QUERY_LIST}"
  exit 2
fi
if ! command -v mafft >/dev/null 2>&1; then
  echo "ERROR: mafft not found on PATH. Install e.g.: conda install -c bioconda mafft"
  exit 3
fi

MOUSE_COL="${MOUSE_COL:-MusMus.pep}"
FOCAL="${FOCAL:-HelKap.pep}"
CLADE3="${CLADE3:-HelKap.pep,HetGla.pep,GeoCap.pep}"

mkdir -p "${OUTDIR}/alignments" "${OUTDIR}/reports"

# Find Orthogroups.tsv
OGTSV=""
if [[ -s "${ORTHO_RESULTS_DIR}/Orthogroups/Orthogroups.tsv" ]]; then
  OGTSV="${ORTHO_RESULTS_DIR}/Orthogroups/Orthogroups.tsv"
else
  # fallback search
  OGTSV="$(find "${ORTHO_RESULTS_DIR}" -maxdepth 3 -type f -name "Orthogroups.tsv" | head -n 1 || true)"
fi
if [[ -z "${OGTSV}" || ! -s "${OGTSV}" ]]; then
  echo "ERROR: could not locate Orthogroups.tsv under: ${ORTHO_RESULTS_DIR}"
  exit 4
fi

# Locate sequence dirs (either/both can exist)
SCO_DIR="${ORTHO_RESULTS_DIR}/Single_Copy_Orthologue_Sequences"
OGSEQ_DIR="${ORTHO_RESULTS_DIR}/Orthogroup_Sequences"

python3 - "${OGTSV}" "${MAPPING_TSV}" "${QUERY_LIST}" "${OUTDIR}" "${MOUSE_COL}" "${SCO_DIR}" "${OGSEQ_DIR}" "${FOCAL}" "${CLADE3}" <<'PY'
import sys, os, re, csv, subprocess
from pathlib import Path
from collections import defaultdict

og_tsv      = Path(sys.argv[1])
map_tsv     = Path(sys.argv[2])
query_list  = Path(sys.argv[3])
outdir      = Path(sys.argv[4])
mouse_col   = sys.argv[5]
sco_dir     = Path(sys.argv[6])
ogseq_dir   = Path(sys.argv[7])
focal       = sys.argv[8]
clade3      = sys.argv[9].split(",")

out_align = outdir / "alignments"
out_rep   = outdir / "reports"
out_align.mkdir(parents=True, exist_ok=True)
out_rep.mkdir(parents=True, exist_ok=True)

def is_prot_id(x: str) -> bool:
    return bool(re.match(r'^[NX]P_\d+(\.\d+)?$', x)) or bool(re.match(r'^[NX]R_\d+(\.\d+)?$', x))

# ---------------------------
# Read mapping: protein_id -> gene
# ---------------------------
prot2gene = {}
gene2prots = defaultdict(set)

with map_tsv.open() as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    # expects columns: protein_id, gene
    for row in rd:
        p = row.get("protein_id","").strip()
        g = row.get("gene","").strip()
        if not p or not g:
            continue
        prot2gene[p] = g
        gene2prots[g].add(p)

# ---------------------------
# Read queries: can be NP/XP or gene symbol
# ---------------------------
wanted_genes = []
wanted_prots = set()
with query_list.open() as fh:
    for line in fh:
        q = line.strip()
        if not q or q.startswith("#"):
            continue
        if is_prot_id(q):
            g = prot2gene.get(q)
            if g:
                wanted_genes.append(g)
                wanted_prots.add(q)
            else:
                # keep as "unknown protein" but still search by raw token in Orthogroups.tsv
                wanted_genes.append(f"UNMAPPED_{q}")
                wanted_prots.add(q)
        else:
            wanted_genes.append(q)
            wanted_prots |= gene2prots.get(q, set())

# unique order-preserving
seen=set()
wanted_genes = [g for g in wanted_genes if not (g in seen or seen.add(g))]

# ---------------------------
# Parse Orthogroups.tsv and find mouse column index
# ---------------------------
with og_tsv.open() as fh:
    header = fh.readline().rstrip("\n").split("\t")
if mouse_col not in header:
    raise SystemExit(f"ERROR: mouse column '{mouse_col}' not found in {og_tsv}. "
                     f"Available columns include: {', '.join(header[:10])} ...")

mouse_idx = header.index(mouse_col)

# map gene -> list of (OG, [mouse_ids_that_match])
gene_hits = defaultdict(list)

# also allow unmapped NP/XP direct match by token search
prot_pattern = re.compile(r'(?:^|, )([A-Za-z0-9_.-]+)')

with og_tsv.open() as fh:
    rd = csv.reader(fh, delimiter="\t")
    hdr = next(rd)
    for row in rd:
        if not row or len(row) <= mouse_idx:
            continue
        og = row[0]
        mouse_cell = row[mouse_idx].strip()
        if not mouse_cell:
            continue
        # OrthoFinder stores gene IDs separated by ", "
        mouse_ids = [x.strip() for x in mouse_cell.split(",") if x.strip()]
        mouse_set = set(mouse_ids)

        # for each wanted gene, check if any of its proteins are in this OG
        for g in wanted_genes:
            if g.startswith("UNMAPPED_"):
                p = g.replace("UNMAPPED_","")
                if p in mouse_set:
                    gene_hits[g].append((og, [p]))
            else:
                prots = gene2prots.get(g, set())
                if not prots:
                    continue
                inter = sorted(prots.intersection(mouse_set))
                if inter:
                    gene_hits[g].append((og, inter))

# write mapping table
map_out = outdir / "gene_orthogroups.tsv"
with map_out.open("w") as oh:
    oh.write("gene\torthogroup\tmouse_ids\n")
    for g in wanted_genes:
        for og, mids in gene_hits.get(g, []):
            oh.write(f"{g}\t{og}\t{','.join(mids)}\n")

# ---------------------------
# Helpers: FASTA IO
# ---------------------------
def read_fasta(p: Path):
    name=None; seq=[]
    with p.open() as fh:
        for line in fh:
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

def write_fasta(p: Path, recs):
    with p.open("w") as oh:
        for n,s in recs:
            oh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                oh.write(s[i:i+60] + "\n")

def find_og_fasta(og: str):
    # prefer single-copy
    p1 = sco_dir / f"{og}.fa"
    if p1.exists() and p1.stat().st_size > 0:
        return p1
    p2 = ogseq_dir / f"{og}.fa"
    if p2.exists() and p2.stat().st_size > 0:
        return p2
    return None

def species_from_id(seqid: str) -> str:
    # your convention: Species_prefix_...
    return seqid.split("_", 1)[0]

# ---------------------------
# For each hit OG: align proteins and scan AA changes
# ---------------------------
summary_rows = []
for gene, hits in gene_hits.items():
    for og, mouse_ids in hits:
        og_fa = find_og_fasta(og)
        if og_fa is None:
            summary_rows.append((gene, og, "NA", "missing_og_fasta", 0, 0))
            continue

        # read OG protein fasta
        recs = list(read_fasta(og_fa))
        if len(recs) < 4:
            summary_rows.append((gene, og, str(og_fa), "too_few_seqs", len(recs), 0))
            continue

        # write a cleaned fasta (first token IDs)
        clean_fa = out_align / f"{gene}__{og}.aa.fa"
        write_fasta(clean_fa, recs)

        aln_fa = out_align / f"{gene}__{og}.aa.aln.fa"
        # run mafft
        cmd = ["mafft", "--auto", str(clean_fa)]
        try:
            mafft_out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL).decode()
        except subprocess.CalledProcessError:
            summary_rows.append((gene, og, str(og_fa), "mafft_failed", len(recs), 0))
            continue
        aln_fa.write_text(mafft_out)

        # parse alignment
        aln = [(n,s) for n,s in read_fasta(aln_fa)]
        seqs = {n:s for n,s in aln}
        spp = {n: species_from_id(n) for n,_ in aln}

        focal_names = [n for n in seqs if spp[n] == focal]
        clade_names = {sp: [n for n in seqs if spp[n] == sp] for sp in clade3}

        # require at least one focal
        if not focal_names:
            summary_rows.append((gene, og, str(og_fa), "no_focal_in_og", len(recs), 0))
            continue
        # for clade3-shared test require all three present (at least one each)
        have_all3 = all(clade_names.get(sp) for sp in clade3)

        # if multiple per species, just take the first (single-copy should be 1 anyway)
        focal_id = focal_names[0]
        clade_id = {sp: clade_names[sp][0] for sp in clade3 if clade_names.get(sp)}

        L = len(seqs[focal_id])
        # basic sanity: all same length
        if any(len(s)!=L for s in seqs.values()):
            summary_rows.append((gene, og, str(og_fa), "aln_len_mismatch", len(recs), 0))
            continue

        # scan positions
        helkap_private = []
        clade3_private = []

        for i in range(L):
            a_f = seqs[focal_id][i]
            if a_f in "-X?":
                continue

            # all other AAs at position (excluding gaps/unknown)
            other = []
            for n,s in seqs.items():
                if n == focal_id:
                    continue
                aa = s[i]
                if aa in "-X?":
                    continue
                other.append(aa)

            # HelKap-only: focal AA not present in any other species
            if other and (a_f not in set(other)):
                helkap_private.append((i+1, a_f, "".join(sorted(set(other)))))

            # clade3-only: HelKap,HetGla,GeoCap share AA, and nobody else has it
            if have_all3:
                aas3 = []
                ok=True
                for sp in clade3:
                    aa = seqs[clade_id[sp]][i]
                    if aa in "-X?":
                        ok=False; break
                    aas3.append(aa)
                if ok and len(set(aas3))==1:
                    shared = aas3[0]
                    # exclude those three and check all remaining
                    other2=[]
                    for n,s in seqs.items():
                        if spp[n] in set(clade3):
                            continue
                        aa=s[i]
                        if aa in "-X?":
                            continue
                        other2.append(aa)
                    if other2 and (shared not in set(other2)):
                        clade3_private.append((i+1, shared, "".join(sorted(set(other2)))))

        # write mutation report
        rep = out_rep / f"{gene}__{og}.mutations.tsv"
        with rep.open("w") as oh:
            oh.write("type\tpos_aa_aln\taa_in_group\taa_in_others\n")
            for pos, aa, oth in helkap_private:
                oh.write(f"{focal}_only\t{pos}\t{aa}\t{oth}\n")
            for pos, aa, oth in clade3_private:
                oh.write(f"{'_'.join(clade3)}_only\t{pos}\t{aa}\t{oth}\n")

        summary_rows.append((gene, og, str(og_fa), "ok", len(recs), len(helkap_private) + len(clade3_private)))

# summary file
sum_out = out_rep / "summary.tsv"
with sum_out.open("w") as oh:
    oh.write("gene\torthogroup\tog_fasta\tstatus\tn_seqs\tn_flagged_sites\n")
    for r in summary_rows:
        oh.write("\t".join(map(str, r)) + "\n")

print("Wrote:", map_out)
print("Alignments:", out_align)
print("Reports:", out_rep)
print("Summary:", sum_out)
PY
