# Selection / OrthoFinder → BPP / PAML / HyPhy pipeline (SLURM)

This repo packages the scripts and conventions from our chats into a
reproducible, versioned workflow.

It's designed for mixed inputs per species:
- genomes (in-house annotation (Mikado) / NCBI-style GFF3 + genome FASTA)
- transcriptomes (Trinity assemblies + curated transcript FASTA)

and produces:
- consistent OrthoFinder-ready proteomes (AA) and matching CDS FASTAs
- OrthoFinder runs + post-run header normalization
- per-OG codon alignments, tested for selection two ways:
  - **PAML** (`codeml`) branch models — one-ratio null vs. two-ratio
    foreground/background LRT
  - **HyPhy** RELAX (relaxed/intensified selection strength) + aBSREL
    (episodic diversifying selection), restricted to the same foreground
    branches as the PAML test so results are directly comparable
- BH-FDR/Holm-corrected, gene-annotated significance tables and summary
  plots for both

## Notice

1) Consistent IDs everywhere

All protein and CDS headers must be prefixed as:

`Species_<original_id>`

This is required because later stages extract CDS by matching OG member IDs
against `ALL.primary.cds.fa`.

2) Avoid shared-file locks on parallel filesystems

The BPP builder writes one status TSV per array task and merges later; do
not use `flock` on shared files.

3) Not every OG has every species

Missing single-copy orthologs is normal. Both the PAML and HyPhy steps
prune the foreground tree down to each OG's *actual* taxon set before
testing it — codeml and HyPhy both require an exact match between tree tips
and alignment sequences, and reusing one static full-taxa tree for every OG
will abort/error on any OG missing even one species. See
`scripts/paml/og_to_codon_paml_codeml.sh` (tree pruning for codeml) and
`scripts/hyphy/prune_hyphy_tree.py` (the same fix for HyPhy).

4) PAML array merges don't wait for other tasks

`scripts/paml/*_array.sh`'s end-of-job summary merge runs "if this is task
1" — which does **not** guarantee task 1 finishes last, so it can produce
an incomplete merge. `scripts/hyphy/merge_hyphy_results.sh` avoids this by
running as a separate job submitted with `--dependency=afterany:<array job
id>` instead — see that script's comments, and consider doing the same for
the PAML merges if you hit this.

## Repository layout

- `scripts/` — genome/transcriptome prep, OrthoFinder helpers, BPP builder
  (unchanged from earlier stages of the pipeline)
  - `scripts/paml/` — codon-alignment building + PAML `codeml` branch-model
    tests (per-baseline build + per-phenotype reuse, significance tables,
    plots)
  - `scripts/hyphy/` — HyPhy RELAX/aBSREL tests on the same alignments and
    foreground branches (tree prep with per-OG pruning, run, SLURM array,
    merge/parse, significance tables, plots)
  - `scripts/annotation/` — genome-wide OG → gene-symbol lookup, built once
    and shared by both `scripts/paml/` and `scripts/hyphy/`
- `slurm/` — sbatch templates for the earlier (OrthoFinder/BPP) stages;
  the PAML/HyPhy stages' own SLURM headers live directly in their
  `scripts/paml/*_array.sh` / `scripts/hyphy/*_array.sh` files instead
- `configs/` — `env.example.sh` (paths/tool locations for `scripts/paml/`
  and `scripts/hyphy/` — copy to `configs/env.sh` and edit) and
  `config.example.yaml` (earlier-stage run parameters)
- `envs/` — conda environment templates
- `workflow/` — optional Snakemake skeleton (thin wrappers around scripts)
- `docs/PIPELINE.md` — every stage, in order, with the script(s) for each
- `legacy/` — superseded early-draft PAML/HyPhy scripts, kept for reference
  only; see `legacy/README.md` for what changed

## Quickstart (minimal, script-driven)

1) Build per-species proteomes + CDS

- Transcriptomes: run `scripts/run_evigene2_trinity_260210.sh <species> <merged_transcripts.fa>`
- Genomes (Mikado): create a primary GFF and extract CDS/proteins with `gffread` using the scripts in `scripts/`

2) Aggregate CDS

Concatenate all per-species `*.cds.fa` files into one:

```bash
cat /path/to/per_species_cds/*.cds.fa > /path/to/run/ALL.primary.cds.fa
```

Generate `imap.txt` (species ↔ individual mapping for BPP):

```bash
bash scripts/make_all_cds_IMap.sh
```

3) Run OrthoFinder

```bash
sbatch scripts/orthofinder.sh
```

4) Normalize OG headers (species prefixes)

```bash
bash scripts/rename_sco_with_species.sh
```

5) Build BPP nucleotide alignments + PHYLIP blocks

```bash
sbatch --export=MIN_TAXA=10,FILTER_STOPS=0,FILTER_FRAME=0,DO_TRIMAL=1 slurm/build_bpp_nt_cds_robust.sh
```

Finish missing OGs: `sbatch slurm/run_missing_ogs.sbatch`, debug with
`scripts/diagnose_missing_ogs.py`.

### 6) Set up config + gene annotation (once per run)

```bash
cp configs/env.example.sh configs/env.sh
# edit configs/env.sh: at minimum set RUN_BASE to your run directory
source configs/env.sh

cd "${SCO_DIR}"   # Single_Copy_Orthologue_Sequences_renamed
bash "${OLDPWD}"/scripts/annotation/extract_og_to_musmus_protein.sh > "${OG2ACC_TSV}"
cd "${OLDPWD}"
```

`MusMus_protein2gene.tsv` (protein accession → gene symbol) is an external
lookup you supply once; point `ACC2GENE_TSV` at it in `configs/env.sh`.

### 7) Build codon alignments + baseline PAML branch model

Full, from-scratch build (mafft + trimAl + codeml) for your baseline
phenotype — needs a PAML `#1`-tagged foreground tree at
`${PAML_DIR}/<phenotype>/species.paml_fg.nwk` first:

```bash
sbatch --export=ALL,PHENOTYPE=only_HelKap scripts/paml/build_codon_alignments_array.sh
```

### 8) PAML branch model for additional phenotypes (reuses step 7's alignments)

```bash
sbatch --export=ALL,PHENOTYPE=acid,TREE_FG=${PAML_DIR}/acid/species.paml_fg.nwk \
    scripts/paml/run_codeml_from_existing_array.sh
```

Then build the significance table, annotate, and plot:

```bash
python3 scripts/paml/build_significance_table.py \
    ${PAML_DIR}/acid/codeml_branch_models_summary.tsv \
    ${PAML_DIR}/acid/acid_significance.tsv

python3 scripts/annotation/annotate_with_gene_names.py \
    ${PAML_DIR}/acid/acid_significance.tsv "${OG2ACC_TSV}" "${ACC2GENE_TSV}" \
    ${PAML_DIR}/acid/acid_significance.annotated.tsv

python3 scripts/paml/plot_paml_summary.py \
    ${PAML_DIR}/acid/acid_significance.annotated.tsv \
    --phenotype acid --correction BH --outdir ${PAML_DIR}/acid
```

### 9) HyPhy RELAX + aBSREL for the same phenotype

```bash
python3 scripts/hyphy/make_hyphy_trees.py \
    ${PAML_DIR}/acid/species.paml_fg.nwk \
    ${HYPHY_DIR}/acid/hyphy_trees/acid

array_jid=$(sbatch --parsable --export=ALL,PHENOTYPE=acid \
    scripts/hyphy/run_hyphy_relax_absrel_array.sh)

sbatch --dependency=afterany:${array_jid} --export=ALL,PHENOTYPE=acid \
    scripts/hyphy/merge_hyphy_results.sh
```

`merge_hyphy_results.sh` runs the full reporting chain (parse → significance
table → gene annotation → plots) in one job — see `docs/PIPELINE.md` for
what each step does.

## Configuration

- `configs/env.example.sh` → copy to `configs/env.sh`, edit `RUN_BASE` (and
  any of the derived paths/tool locations you need to override), then
  `source configs/env.sh` before running anything in `scripts/paml/` or
  `scripts/hyphy/`. Every script also accepts the same variables as
  per-invocation overrides (`RUN_BASE=... sbatch --export=ALL ...`).
- `configs/config.example.yaml` — parameters for the earlier
  (transcriptome/genome/OrthoFinder/BPP) stages.

Scripts in `scripts/paml/` and `scripts/hyphy/` locate their own sibling
scripts (an array driver finding its per-OG worker, `merge_hyphy_results.sh`
finding `scripts/annotation/`) relative to their own location on disk, not
a hardcoded path — so each `scripts/<stage>/` directory is self-contained
and works wherever you clone or copy it, including copying individual
files out to a flat `~/scripts/`-style directory on a cluster.
