# Pipeline stages

1) Transcriptomes → primary_by_gene (Evigene + TransDecoder)
   - `scripts/run_evigene2_trinity_260210.sh`

2) Genomes (Mikado/NCBI GFF) → primary transcripts
   - `scripts/mikado_primary_gff.py` or `scripts/keep_one_tx_per_gene.py`
   - `gffread` extraction + header prefixing (`scripts/make_cds_fromGFF_*.sh`)

3) OrthoFinder run
   - `scripts/orthofinder.sh`

4) Post-OrthoFinder header normalization
   - `scripts/rename_sco_with_species.sh`

5) OG → nucleotide alignment → BPP PHYLIP
   - `scripts/og_to_codon2_phy.sh`
   - `scripts/build_bpp_nt_cds_robust.sh` (array driver)
   - `scripts/diagnose_missing_ogs.py` + `slurm/run_missing_ogs.sbatch`

6) Genome-wide gene-symbol annotation (build once, reuse everywhere)
   - `scripts/annotation/extract_og_to_musmus_protein.sh` — every
     single-copy OG's MusMus protein accession → `OG_to_MusMus_protein.tsv`
   - `scripts/annotation/annotate_with_gene_names.py` — joins that against
     `MusMus_protein2gene.tsv` to add a `gene_symbol` column to any
     OG-keyed table. Run once per OrthoFinder run; every PAML and HyPhy
     phenotype/test below reuses the same two lookup files.

7) OG → codon alignments for PAML + HyPhy (shared step)
   - `scripts/paml/og_to_codon_paml_codeml.sh` — per-OG worker: filters CDS
     (frame, internal stops), aligns (mafft), optionally trims (trimAl),
     backtranslates to a codon alignment, writes both PAML PHYLIP
     (`paml_phy/`) and HyPhy FASTA (`hyphy_fa/`, byte-identical content) for
     that OG, and (if `RUN_CODEML=1`) prunes the foreground tree to that
     OG's actual taxon set and runs codeml immediately.
   - `scripts/paml/build_codon_alignments_array.sh` — SLURM array driver:
     runs the worker above over every single-copy OG for one **baseline**
     phenotype (default `only_HelKap`). This is the expensive, from-scratch
     step (mafft + trimAl + codeml) — run it once per baseline, then reuse
     its alignments for every other phenotype (step 8) rather than
     rebuilding them.

8) PAML branch-model tests, per phenotype
   - `scripts/paml/codeml_from_existing_alignment.sh` — per-OG worker:
     reuses an existing `paml_phy/<OG>.codon.phy` from step 7, re-prunes the
     tree to that OG's taxa against a **new** phenotype's foreground tree,
     and reruns just codeml (no realignment).
   - `scripts/paml/run_codeml_from_existing_array.sh` — SLURM array driver
     for the worker above. Run once per additional phenotype
     (`acid`, `acid_wBat`, `capsaicin`, `HelKap_HetGla`, ...).
   - `scripts/paml/build_significance_table.py` — BH-FDR + Holm-Bonferroni
     correction over a phenotype's merged codeml summary.
   - `scripts/annotation/annotate_with_gene_names.py` — adds gene symbols
     (step 6's lookup) to the significance table.
   - `scripts/paml/plot_paml_summary.py` — p-value histogram + Manhattan-
     style scan, top hits labeled by gene name.

9) HyPhy RELAX + aBSREL tests, per phenotype
   Reuses the SAME `hyphy_fa/<OG>.codon.fasta` files from step 7 (no
   realignment) and tests the SAME Test/foreground branches as that
   phenotype's PAML tree, so results are directly comparable to step 8's
   codeml LRT.
   - `scripts/hyphy/make_hyphy_trees.py` — converts a phenotype's PAML
     `#1`-tagged foreground tree into HyPhy-labeled trees: every leaf
     `{Test}`/`{Reference}` (for RELAX) or just the foreground leaves
     `{Test}` (for aBSREL, restricting the scan to `--branches Test`
     instead of a genome-wide all-branch scan).
   - `scripts/hyphy/prune_hyphy_tree.py` — prunes a HyPhy-labeled tree down
     to one OG's actual taxon set (same problem/fix as the codeml tree
     pruning in step 7 — not every OG has every species, and HyPhy requires
     an exact match between tree tips and alignment sequences).
   - `scripts/hyphy/run_hyphy_relax_absrel.sh` — per-OG-list worker: for
     each OG, prunes both trees (via the script above) and runs
     `hyphy relax` and/or `hyphy absrel` (`ONLY=relax|absrel|both`).
   - `scripts/hyphy/run_hyphy_relax_absrel_array.sh` — SLURM array driver:
     splits the phenotype's full OG list into interleaved chunks and hands
     each to the worker above. Genome-wide by default (not just a
     significant subset).
   - `scripts/hyphy/merge_hyphy_results.sh` — run as a **separate job**
     with `--dependency=afterany:<array job id>` after the array above
     finishes (see the script's own comments for why this isn't folded
     into "task 1" the way the PAML array merge is). Runs the full
     reporting chain in one go:
     `parse_hyphy_results.py` → `scripts/hyphy/build_hyphy_significance_table.py`
     → `scripts/annotation/annotate_with_gene_names.py` →
     `scripts/hyphy/plot_hyphy_summary.py` (run twice: `--test relax`,
     `--test absrel` — these are two different hypothesis tests, corrected
     and plotted independently).

10) Mouse-target gene mapping and mutation scan
    - `scripts/extract_mouse_geneogs_scan.py`
    - `scripts/extract_mouse_gene_ogs_and_scan.sh`

See `legacy/README.md` for what changed between the old `codeml_batch_launcher.sh`
/ `og_to_codon_paml_hyphy.sh` / `run_absrel_fg.sbatch` and the current
`scripts/paml/` + `scripts/hyphy/` pipeline (steps 7-9 above).
