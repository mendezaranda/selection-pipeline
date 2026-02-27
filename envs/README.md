Put your conda env YAMLs here.

Recommended split:
- `orthofinder.yaml`: orthofinder + diamond
- `trinity.yaml`: trinity + transdecoder + cd-hit + perl deps for evigene
- `phylo.yaml`: mafft + trimal + hyphy + paml (codeml)

Cluster modules differ; keep sbatch templates in `slurm/`.
