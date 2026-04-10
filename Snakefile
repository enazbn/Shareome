configfile: "config/config.yaml"

import os

FASTAS = config["fastas"]
INPUT_DIR = config["input_dir"]
OUTPUT_DIR = config["output_dir"]

HUMAN_NAME = "human_refseq_nr_preclean1"

VIRUS_NAMES = [
    os.path.splitext(f)[0]
    for f in FASTAS
    if os.path.splitext(f)[0] != HUMAN_NAME
]

ALL_NAMES = [os.path.splitext(f)[0] for f in FASTAS]

rule all:
    input:
        expand(OUTPUT_DIR + "/{name}_kmers.csv", name=ALL_NAMES),
        expand(OUTPUT_DIR + "/{name}_matched.csv", name=VIRUS_NAMES),
        expand(OUTPUT_DIR + "/{name}_matched.parquet", name=VIRUS_NAMES),
        expand(OUTPUT_DIR + "/{name}_peptide_core.parquet", name=VIRUS_NAMES)
rule kmerslice:
    input:
        fasta=lambda wc: INPUT_DIR + "/" + wc.name + ".fasta"
    output:
        csv=OUTPUT_DIR + "/{name}_kmers.csv"
    log:
        "logs/{name}_kmers.log"
    conda:
        "workflow/envs/kmerslicer.yaml"
    shell:
        """
        python tools/kmerslicer/kmerslicer.py \
        {input.fasta} \
        {output.csv} \
        -k {config[k]} \
        --format csv \
        --skip-ambiguous \
        > {log} 2>&1
        """

rule match_kmers:
    input:
        human=OUTPUT_DIR + "/" + HUMAN_NAME + "_kmers.csv",
        virus=OUTPUT_DIR + "/{name}_kmers.csv"
    output:
        matched=OUTPUT_DIR + "/{name}_matched.csv"
    log:
        "logs/{name}_matched.log"
    conda:
        "workflow/envs/kmerslicer.yaml"
    shell:
        """
        python workflow/scripts/match_kmers.py \
        {input.human} \
        {input.virus} \
        {output.matched} \
        > {log} 2>&1
        """
rule csv_to_parquet:
    input:
        csv=OUTPUT_DIR + "/{name}_matched.csv"
    output:
        parquet=OUTPUT_DIR + "/{name}_matched.parquet"
    log:
        "logs/{name}_parquet.log"
    conda:
        "workflow/envs/kmerslicer.yaml"
    shell:
        """
        python workflow/scripts/csv_to_parquet.py \
        {input.csv} \
        {output.parquet} \
        > {log} 2>&1
        """
rule aggregate_peptide_core:
    input:
        parquet=OUTPUT_DIR + "/{name}_matched.parquet"
    output:
        core=OUTPUT_DIR + "/{name}_peptide_core.parquet"
    log:
        "logs/{name}_peptide_core.log"
    conda:
        "workflow/envs/kmerslicer.yaml"
    shell:
        """
        python workflow/scripts/aggregate_matched_parquet.py \
        {input.parquet} \
        {output.core} \
        > {log} 2>&1
        """
