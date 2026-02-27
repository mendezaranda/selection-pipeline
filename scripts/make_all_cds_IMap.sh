CDS_DIR="/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/02_cds"
OUTDIR="/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/08_bpp"

mkdir -p "$OUTDIR"
#cat "${CDS_DIR}"/*.cds.fa > "${OUTDIR}/ALL.primary.cds.fa"

# Make imap.txt: individualID speciesName (individualID is after ^ in alignments, here Species)
python3 - <<PY
import glob
cds_dir="${CDS_DIR}"
species=set()
for fn in glob.glob(cds_dir+"/*.cds.fa"):
    sp=fn.split("/")[-1].split(".")[0]  # Species from Species.primary...
    sp=sp.split("_")[0]                 # safe if any underscores appear
    species.add(sp)
with open("${OUTDIR}/imap.txt","w") as o:
    for sp in sorted(species):
        o.write(f"{sp}\t{sp}\n")
print("Wrote imap.txt with", len(species), "species")
PY
