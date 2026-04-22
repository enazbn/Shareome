configfile: "config/config.yaml"

import os

FASTAS = config["fastas"]
INPUT_DIR = config["input_dir"]
OUTPUT_DIR = config["output_dir"]
K = config["k"]
MISMATCH_THRESHOLDS = config["mismatch_thresholds"]

HUMAN_NAME = "human_refseq_nr_cleaned"

ALL_NAMES = [os.path.splitext(f)[0] for f in FASTAS]

VIRUS_NAMES = [
    os.path.splitext(f)[0]
    for f in FASTAS
    if os.path.splitext(f)[0] != HUMAN_NAME
]


rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "{name}_kmers.csv"), name=ALL_NAMES),
        expand(
            os.path.join(OUTPUT_DIR, "{name}_mm{mm}_matched.csv"),
            name=VIRUS_NAMES,
            mm=MISMATCH_THRESHOLDS
        ),
        expand(
            os.path.join(OUTPUT_DIR, "{name}_mm{mm}_matched.parquet"),
            name=VIRUS_NAMES,
            mm=MISMATCH_THRESHOLDS
        ),
        expand(
            os.path.join(OUTPUT_DIR, "{name}_mm{mm}_peptide_core.parquet"),
            name=VIRUS_NAMES,
            mm=MISMATCH_THRESHOLDS
        )


rule kmerslice:
    input:
        fasta=lambda wc: os.path.join(INPUT_DIR, wc.name + ".fasta")
    output:
        csv=os.path.join(OUTPUT_DIR, "{name}_kmers.csv")
    log:
        "logs/{name}_kmers.log"
    conda:
        "workflow/envs/kmerslicer.yaml"
    shell:
        """
        python tools/kmerslicer/kmerslicer.py \
            {input.fasta} \
            {output.csv} \
            -k {K} \
            --format csv \
            --skip-ambiguous \
            > {log} 2>&1
        """


rule match_kmers_threshold:
    input:
        human=os.path.join(OUTPUT_DIR, HUMAN_NAME + "_kmers.csv"),
        virus=os.path.join(OUTPUT_DIR, "{name}_kmers.csv")
    output:
        matched=os.path.join(OUTPUT_DIR, "{name}_mm{mm}_matched.csv")
    params:
        mm=lambda wc: int(wc.mm)
    log:
        "logs/{name}_mm{mm}_matched.log"
    conda:
        "workflow/envs/kmerslicer.yaml"
    shell:
        """
        python workflow/scripts/match_kmers_threshold.py \
            {input.human} \
            {input.virus} \
            {params.mm} \
            {output.matched} \
            > {log} 2>&1
        """


rule csv_to_parquet:
    input:
        csv=os.path.join(OUTPUT_DIR, "{name}_mm{mm}_matched.csv")
    output:
        parquet=os.path.join(OUTPUT_DIR, "{name}_mm{mm}_matched.parquet")
    log:
        "logs/{name}_mm{mm}_parquet.log"
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
        parquet=os.path.join(OUTPUT_DIR, "{name}_mm{mm}_matched.parquet")
    output:
        core=os.path.join(OUTPUT_DIR, "{name}_mm{mm}_peptide_core.parquet")
    log:
        "logs/{name}_mm{mm}_peptide_core.log"
    conda:
        "workflow/envs/kmerslicer.yaml"
    shell:
        """
        python workflow/scripts/aggregate_matched_parquet.py \
            {input.parquet} \
            {output.core} \
            > {log} 2>&1
        """