#!/bin/bash

# mikado_make_cds.sh
set -eo pipefail
species="${1:?species code}"
gff="${2:?mikado gff3}"
genome="${3:?genome fasta}"
out=/fast/AG_Lewin/dmendez/selection_pipeline/2026-02-10_run01/02_cds/${species}_cds.fa

tmp="$(mktemp -p . ${species}.XXXX.cds.fa)"
gffread -g "$genome" -x "$tmp" "$gff"

awk -v sp="$species" '
  /^>/{
    id=substr($0,2); sub(/[[:space:]].*$/, "", id);
    print ">" sp "_" id; next
  }
  {print}
' "$tmp" > "$out"
rm -f "$tmp"

## usage: bash ~/scripts/make_cds_fromGFF_260210.sh BatSui /fast/AG_Lewin/dmendez/genomes/BatSui/BatSui_Dustin/BatSui.edited.geneSymbol.gff /fast/AG_Lewin/dmendez/genomes/BatSui/BatSui_Dustin/BatSui.pre.a4bm0l.bp.p_ctg.cleaned.sm.fa 02_cds/BatSui.cds.fasta