import os
import glob

POS,FRR = glob_wildcards("raw_reads/{pos}_{frr}_001.fastq.gz")

# Véritable output voulu
rule all:
    input:
        expand("fastqc/{pos}_{frr}_001_fastqc.{extension}", pos=POS, frr=FRR, extension=["zip", "html"]),
        expand("assemblyqc/{pos}_assemblyqc/report.html", pos=POS),
        expand("genome_annotation/{pos}_annotation/{pos}.tab", pos=POS)

# FastQC
rule fastQC:
    input:
        rawreads = "raw_reads/{pos}_{frr}_001.fastq.gz"
    output:
        zip = "fastqc/{pos}_{frr}_001_fastqc.zip",
        html = "fastqc/{pos}_{frr}_001_fastqc.html"
    message:
        "Make QC of raw reads."
    threads: 2
    params:
        out_prefix = "fastqc/"
    shell:
        """
        fastqc -o {params.out_prefix} -t {threads} {input.rawreads}
        """

# Trimmed reads
rule trimmomatic:
    input:
        r1 = "raw_reads/{pos}_R1_001.fastq.gz",
        r2 = "raw_reads/{pos}_R2_001.fastq.gz",
    output:
        forward_paired = "trimmed_reads/{pos}_1P.fastq",
        reverse_paired = "trimmed_reads/{pos}_2P.fastq",
        forward_unpaired = "trimmed_reads/{pos}_1U.fastq",
        reverse_unpaired = "trimmed_reads/{pos}_2U.fastq"
    message:
        "Trimm reads."
    threads: 8
    log:
        "trimmed_reads/{pos}.log"
    params:
        basename = "trimmed_reads/{pos}.fastq"
    shell:
        """
        trimmomatic PE -threads {threads} \
        {input.r1} {input.r2} \
        -baseout {params.basename} \
        ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:8:keepBothReads \
        HEADCROP:15 SLIDINGWINDOW:4:25 \
        2>{log}
        """

# Mix unpaired read
rule unpaired:
    input:
        forward_unpaired = rules.trimmomatic.output.forward_unpaired,
        reverse_unpaired = rules.trimmomatic.output.reverse_unpaired
    output:
        forward_and_reverse_unpaired = "trimmed_reads/{pos}_unpaired.fastq"
    message:
        "Mix unpaired read."
    shell:
        """
        cat {input.forward_unpaired} {input.reverse_unpaired} > {output.forward_and_reverse_unpaired}
        """

# Assembly genome
rule assembly:
    input:
        r1 = rules.trimmomatic.output.forward_paired,
        r2 = rules.trimmomatic.output.reverse_paired,
        unpaired = rules.unpaired.output.forward_and_reverse_unpaired
    output:
        genome_assembly = "genome_assembly/{pos}_assembly/assembly.fasta"
    message:
        "Assembly genome."
    log:
        "genome_assembly/{pos}_assembly/unicycler.log"
    threads: 32
    params:
        out_prefix = "genome_assembly/{pos}_assembly/"
    shell:
        """
        unicycler -1 {input.r1} \
        -2 {input.r2} \
        -s {input.unpaired} \
        -o {params.out_prefix} \
        --verbosity 2 \
        --min_fasta_length 200 \
        2>{log} \
        -t {threads}
        """    

# Genome annotation
rule prokka:
    input:
        genome_assembly = rules.assembly.output.genome_assembly
    output:
        genome_annotation = "genome_annotation/{pos}_annotation/{pos}.gff"
    message:
        "Annotate the genome."
    log:
        "genome_annotation/{pos}_annotation/prokka.log"
    params:
        out_prefix = "genome_annotation/{pos}_annotation/",
        filename_prefix = "{pos}"
    threads: 32
    shell:
        """
        prokka --outdir {params.out_prefix} \
        --prefix {params.filename_prefix} \
        --genus 'Escherichia' \
        --species 'coli' \
        --cpus {threads} \
        --force \
        2>{log} \
        {input.genome_assembly}
        """

# Assembly QC
rule quast:
    input:
        r1 = rules.trimmomatic.output.forward_paired,
        r2 = rules.trimmomatic.output.reverse_paired,
        unpaired = rules.unpaired.output.forward_and_reverse_unpaired,
        genome_annotation = rules.prokka.output.genome_annotation,
        genome_assembly = rules.assembly.output.genome_assembly

    output:
        html = "assemblyqc/{pos}_assemblyqc/report.html"
    message:
        "Assembly QC"
    log:
        "assemblyqc/{pos}_assemblyqc/quast.log"
    params:
        out_prefix = "assemblyqc/{pos}_assemblyqc/"
    threads: 32
    shell:
        """
        quast -o {params.out_prefix} \
        -g {input.genome_assembly} \
        -t {threads} \
        -1 {input.r1} \
        -2 {input.r2} \
        --single {input.unpaired} \
        --gene-thresholds 0,1000 \
        {input.genome_assembly} \
        2>{log} \
        --glimmer
        """

# Research for antimicrobial genes resistance
rule abricate:
    input:
        genome_annotation = rules.prokka.output.genome_annotation
    output:
        tab_amrg = "genome_annotation/{pos}_annotation/{pos}.tab"
    message:
        "Research for antimicrobial genes resistance"
    threads: 32
    shell:
        """
        abricate --threads {threads} {input.genome_annotation} > {output.tab_amrg}
        """