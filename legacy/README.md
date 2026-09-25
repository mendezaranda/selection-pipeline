# Legacy scripts

These are the original, early-draft versions of the PAML/HyPhy steps,
superseded by `scripts/paml/` and `scripts/hyphy/`. Kept for reference only
-- do not use for a new run.

What changed since these were written (see the main README and
`docs/PIPELINE.md` for the current pipeline):

- **Per-OG tree pruning.** Not every orthogroup has all species (missing
  single-copy orthologs is normal). `codeml_batch_launcher.sh` and
  `og_to_codon_paml_hyphy.sh` reuse one static tree for every OG, which
  aborts on any OG missing a taxon. `scripts/paml/og_to_codon_paml_codeml.sh`
  and `scripts/hyphy/prune_hyphy_tree.py` fix this by pruning the tree to
  each OG's actual taxon set before running codeml/HyPhy.
- **RELAX added alongside aBSREL**, both restricted to the phenotype's own
  foreground branches (`--branches Test` / `--test Test --reference
  Reference`) instead of `run_absrel_fg.sbatch`'s plain all-branches
  `--tree` call.
- **Incremental, resumable runs.** The current scripts skip any OG whose
  output already exists (a parseable codeml lnL, or a non-empty HyPhy JSON),
  so a rerun after a bugfix or a timeout only reprocesses what didn't
  actually finish -- not a full restart.
- **Alignment reuse across phenotypes.** `codeml_from_existing_alignment.sh`
  reruns only the tree-pruning + codeml step against a new phenotype's
  foreground tree, reusing the (expensive) mafft/trimAl alignment from a
  prior full run instead of rebuilding it from scratch per phenotype.
- **Significance tables + plots.** `build_significance_table.py` /
  `build_hyphy_significance_table.py` add BH-FDR and Holm-Bonferroni
  corrections; `plot_paml_summary.py` / `plot_hyphy_summary.py` add
  diagnostic p-value histograms and Manhattan-style scans with genome-wide
  gene-symbol annotation (`scripts/annotation/`), none of which existed in
  this legacy version.
