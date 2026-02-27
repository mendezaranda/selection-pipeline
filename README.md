# Selection / OrthoFinder → BPP / PAML / HyPhy pipeline (SLURM)

This repo packages the scripts and conventions from our chat into a reproducible, versioned workflow.

It’s designed for mixed inputs per species:
- genomes (Mikado / NCBI style GFF3 + genome FASTA)
- transcriptomes (Trinity assemblies + curated transcript FASTA)

and produces:
- consistent OrthoFinder-ready proteomes (AA) and matching CDS FASTAs
- OrthoFinder runs + post-run header normalization
- per-OG nucleotide and codon alignments suitable for BPP / ASTRAL / PAML / HyPhy (aBSREL)

## Non-negotiable conventions

1) Consistent IDs everywhere

All protein and CDS headers must be prefixed as:

`Species_<original_id>`

This is required because later stages extract CDS by matching OG member IDs against `ALL.primary.cds.fa`.

2) Avoid fragile interactive shell code in batch jobs

If your `~/.bashrc` prints quota or runs commands that can fail under `set -e`, guard it:

```bash
# in ~/.bashrc
if [[ $- == *i* ]]; then
  # interactive-only stuff
fi
```

3) Avoid shared-file locks on parallel filesystems

The BPP builder writes one status TSV per array task and merges later; do not use `flock` on shared files.

## Repository layout

- `scripts/` core per-step scripts (genome, transcriptome, OrthoFinder helpers, OG → alignments)
- `slurm/` sbatch templates for cluster execution
- `configs/` example configs (paths, species lists, parameters)
- `envs/` conda environment templates
- `workflow/` optional Snakemake skeleton (thin wrappers around scripts)

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

Use the provided sbatch template as a starting point:

```bash
sbatch scripts/orthofinder.sh
```

4) Normalize OG headers (species prefixes)

```bash
bash scripts/rename_sco_with_species.sh
```

5) Build BPP nucleotide alignments + PHYLIP blocks

Option A (stride array over all OGs):

```bash
sbatch --export=MIN_TAXA=10,FILTER_STOPS=0,FILTER_FRAME=0,DO_TRIMAL=1 slurm/build_bpp_nt_cds_robust.sh
```

Option B (finish missing OGs):

```bash
sbatch slurm/run_missing_ogs.sbatch
```

Debug missing OGs:

```bash
python3 scripts/diagnose_missing_ogs.py \
  --sco_dir /path/to/Single_Copy_Orthologue_Sequences_renamed \
  --all_cds /path/to/ALL.primary.cds.fa \
  --missing_list missing_ogs.txt \
  --out_tsv diagnose.tsv \
  --min_taxa 10
```

6) Build codon alignments for PAML/HyPhy (optional)

```bash
# per OG
bash scripts/og_to_codon_paml_hyphy.sh OG0000001.fa ALL.primary.cds.fa outdir 16
```

7) Run codeml / aBSREL (optional)

- codeml batch launcher: `scripts/codeml_batch_launcher.sh`
- HyPhy aBSREL sbatch template: `slurm/run_absrel_fg.sbatch`

## Configuration

Start from `configs/config.example.yaml` and adapt paths to your run.

## Notes on hardcoded paths

Some scripts currently include lab-specific absolute paths (e.g. `/fast/...`, `/home/dmendez/...`).
Before publishing or sharing broadly:
- replace these with config variables or CLI arguments
- keep sbatch templates in `slurm/` as examples, not required defaults

## License

Add a LICENSE before publishing publicly.
