Shareome Snakemake Pipeline

A modular and scalable Snakemake workflow for protein k-mer generation and cross-proteome matching.

This pipeline is part of the **Shareome framework**, aimed at systematically identifying shared peptide sequences between viral and human proteomes.

This workflow:

1. Takes multiple protein FASTA files (human + viruses)
2. Generates position-resolved k-mers using a custom k-mer slicer
3. Produces structured CSV outputs for downstream matching and analysis

The pipeline is designed to be:
- Modular (step-by-step execution)
- Reproducible (Conda environments)
- Scalable (laptop → HPC)
- Analysis-ready (CSV → Parquet workflows)

