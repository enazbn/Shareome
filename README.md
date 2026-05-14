# Shareome Snakemake Pipeline

A modular and scalable Snakemake workflow for protein k-mer generation and cross-proteome peptide similarity analysis.

This pipeline is part of the Shareome framework, which aims to systematically identify short shared or similar peptide sequences between viral and human proteomes. The workflow supports both identity-based matching and substitution-matrix-based similarity scoring, allowing users to explore exact matches, near matches, and biochemically similar peptide relationships.

---

## Overview

The Shareome pipeline:

1. Takes multiple protein FASTA files, including one human reference proteome and one or more viral proteomes
2. Generates position-resolved protein k-mers using a custom k-mer slicer
3. Performs cross-proteome matching between viral and human k-mers
4. Supports two mutually exclusive matching modes:
   - **Mismatch mode**: exact and near-exact matching based on the number of amino-acid mismatches
   - **BLOSUM mode**: similarity-based matching using a substitution matrix such as BLOSUM62
5. Converts large CSV outputs to Parquet format
6. Aggregates match results into peptide-level core tables for downstream analysis

The pipeline is designed to be:

- **Modular**: each step can be run, inspected, and modified independently
- **Reproducible**: compatible with Conda-based Snakemake environments
- **Scalable**: suitable for laptop-scale testing and HPC-scale execution
- **Analysis-ready**: produces structured CSV and Parquet outputs for downstream filtering, scoring, and visualization

---

## Pipeline structure

```text
FASTA files
   ↓
k-mer slicing
   ↓
position-resolved k-mer CSV files
   ↓
cross-proteome matching
   ↓
matched peptide CSV files
   ↓
CSV to Parquet conversion
   ↓
peptide-level aggregation
   ↓
analysis-ready peptide core tables