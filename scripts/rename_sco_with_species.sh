#!/bin/bash

set -eo pipefail

# Usage:
#   bash rename_sco_with_species.sh /path/to/03b_orthofinder/OrthoFinder/Results_XXXX
#
# It will create:
#   <RESULTS_DIR>/Single_Copy_Orthologue_Sequences_renamed

#RESULTS_DIR="${1:?Provide OrthoFinder Results directory (contains Single_Copy_Orthologue_Sequences)}"
RESULTS_DIR="/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/03b_orthofinder/out/Results_Feb12"
SCO_DIR="${RESULTS_DIR}/Orthogroup_Sequences"
OUT_DIR="${RESULTS_DIR}/Orthogroup_Sequences_renamed"
LOG="${OUT_DIR}/rename_ogs.log"

mkdir -p "${OUT_DIR}"
: > "${LOG}"

if [[ ! -d "${SCO_DIR}" ]]; then
  echo "ERROR: missing ${SCO_DIR}" >&2
  exit 2
fi

python3 - "$SCO_DIR" "$OUT_DIR" "$LOG" <<'PY'
import re
import sys
from pathlib import Path

sco_dir = Path(sys.argv[1])
out_dir = Path(sys.argv[2])
log_fn  = Path(sys.argv[3])

re_batsui = re.compile(r'^mikado(?:$|[_\.\-])')
re_fukans = re.compile(r'^g\d+\.t\d+$')
re_prefixed = re.compile(r'^[A-Za-z][A-Za-z0-9]{2,}_.+')

def fasta_iter(p: Path):
    h = None
    s = []
    with p.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(s)
                h = line[1:].split()[0]  # first token only
                s = []
            else:
                s.append(line.strip())
        if h is not None:
            yield h, "".join(s)

def write_fasta(records, outp: Path):
    with outp.open("w") as oh:
        for h, seq in records:
            oh.write(f">{h}\n")
            for i in range(0, len(seq), 60):
                oh.write(seq[i:i+60] + "\n")

og_files = sorted(list(sco_dir.glob("OG*.fa")) + list(sco_dir.glob("OG*.fasta")))
if not og_files:
    raise SystemExit(f"ERROR: no OG*.fa or OG*.fasta found in {sco_dir}")

with log_fn.open("w") as lg:
    lg.write(f"SCO input:\t{sco_dir}\n")
    lg.write(f"SCO output:\t{out_dir}\n")
    lg.write("file\tseqs\trenamed_BatSui\trenamed_FukAns\tkept_prefixed\tkept_other\n")

for fp in og_files:
    n = n_b = n_f = n_pref = n_other = 0
    out_records = []
    for h, seq in fasta_iter(fp):
        n += 1
        if re_prefixed.match(h):
            n_pref += 1
            newh = h
        elif re_batsui.match(h):
            n_b += 1
            newh = "BatSui_" + h
        elif re_fukans.match(h):
            n_f += 1
            newh = "FukAns_" + h
        else:
            n_other += 1
            newh = h
        out_records.append((newh, seq))

    outp = out_dir / fp.name
    write_fasta(out_records, outp)

    with log_fn.open("a") as lg:
        lg.write(f"{fp.name}\t{n}\t{n_b}\t{n_f}\t{n_pref}\t{n_other}\n")

print("Wrote renamed SCO OG fastas to:", out_dir)
print("Log:", log_fn)
PY

echo "Done."
echo "Output: ${OUT_DIR}"
echo "Log: ${LOG}"
