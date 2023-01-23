import os
import glob

POS,FRR = glob_wildcards("raw_reads/{pos}_{frr}_001.fastq.gz")

folder = ["fastqc", "trimmed_reads", "snippy"]
for f in folder:
    os.mkdir(f)

# Rule for the latest files wanted
rule all:
    input:
        expand("fastqc/{pos}_{frr}_001_fastqc.{extension}", pos=POS, frr=FRR, extension=["zip", "html"]),
        "snippy/output/variants.tab"

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
        r2 = "raw_reads/{pos}_R2_001.fastq.gz"
    output:
        forward_paired = "trimmed_reads/{pos}_1P.fastq",
        reverse_paired = "trimmed_reads/{pos}_2P.fastq"
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

# Variant calling need genome reference
rule snippy:
    input:
        forward_paired = rules.trimmomatic.output.forward_paired,
        reverse_paired = rules.trimmomatic.output.reverse_paired,
        genome_ref = config["reference"]
    output:
        snps_tab = "snippy/{pos}/{pos}_snps.tab"
    message:
        "Call the variant"
    threads: 32
    log:
        "snippy/{pos}/snippy.log"
    params:
        outdir = "snippy/{pos}/",
        prefix = "{pos}_snps"
    shell:
        """
        snippy --cpus {threads} \
        --force \
        --outdir {params.outdir} \
        --prefix {params.prefix} \
        --ref {input.genome_ref} \
        --R1 {input.forward_paired} \
        --R2 {input.reverse_paired} \
        2>{log}
        """

rule take_tab_path:
    input:
        snps_tab = expand("snippy/{pos}/{pos}_snps.tab", pos=POS)
    output:
        snps_path_tab = "snippy/snps_path_tab.txt"
    message:
        "Regroup all the path of the snps.tab"
    shell:
        """
        realpath {input.snps_tab} >> {output.snps_path_tab}
        """

rule regroup_variant:
    input:
        snps_path_tab = rules.take_tab_path.output.snps_path_tab,
        metadata = "metadata.csv"
    output:
        all_variants = "snippy/output/variants.tab"
    message:
        "Regroup all the variants according to the metadata"
    params:
        directory = "snippy/"
    shell:
        """ 
        python3 snp_regroup.py {input.metadata} {input.snps_path_tab} {params.directory}
        """
