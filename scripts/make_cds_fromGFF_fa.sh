species=HelKap
gff=/fast/AG_Lewin/dmendez/genomes/HelPot/Kapiti_annotations_Oct2025/full_annotation.edited.geneSymbol.gff
genome=/fast/AG_Lewin/dmendez/genomes/HelPot/earlgrey2/Heliophobius_argenteocinereus_EarlGrey/Heliophobius_argenteocinereus_summaryFiles/Heliophobius_argenteocinereus.softmasked.fasta
tmp="/fast/AG_Lewin/dmendez/orthofinder/20251219/primary/input_clean/cds/${species}.tmp.cds.fa"
out="/fast/AG_Lewin/dmendez/orthofinder/20251219/primary/input_clean/cds/${species}.primary.cds.fa"

# Extract all CDS sequences per transcript ID (mRNA)
gffread -g "$genome" -x "$tmp" "$gff"

# Prefix to match OrthoFinder protein IDs:
awk -v sp="$species" '
  /^>/{sub(/^>/, ">" sp "_", $0); print; next}
  {print}
' "$tmp" > "$out"


# species=BatSui
# OG=OG0011139
# grep '^>' /fast/AG_Lewin/dmendez/orthofinder/20251219/primary/input_clean/OrthoFinder/Results_Jan13/Single_Copy_Orthologue_Sequences/$OG.fa | sed 's/^>//' | awk '{print $1}' | grep "^${species}_" | sort > ${species}_og.ids
# grep '^>' /fast/AG_Lewin/dmendez/orthofinder/20251219/primary/input_clean/cds/${species}.primary.cds.fa | sed 's/^>//' | awk '{print $1}' | sort > ${species}_cds.ids
# comm -23 ${species}_og.ids ${species}_cds.ids | head
