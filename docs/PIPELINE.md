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

6) OG → codon alignments for PAML/HyPhy
   - `scripts/og_to_codon_paml_hyphy.sh`

7) Downstream selection tests
   - `scripts/codeml_batch_launcher.sh`
   - `slurm/run_absrel_fg.sbatch`

8) Mouse-target gene mapping and mutation scan
   - `scripts/extract_mouse_geneogs_scan.py`
   - `scripts/extract_mouse_gene_ogs_and_scan.sh`
