import os
from snakemake.utils import min_version
min_version("6.2")


report: "report/workflow.rst"


FASTQ=["SRR13435231", "SRR13435229"]


localrules: all, sample_sheet, trimmed_reads, unzip
rule all:
    input: expand(["reads/{accession}_1.fastq.gz", "reads/{accession}_2.fastq.gz", "results/Assembly/MEGAHIT/{accession}.contigs.fa", "results/metator/{accession}/bin_summary.txt", "results/virsorter2/{accession}/final-viral-score.tsv", "results/metator/{accession}/contact_map"], accession=FASTQ)


rule get_fastq_pe_gz:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        "reads/{accession}_1.fastq.gz",
        "reads/{accession}_2.fastq.gz",
    log:
        "logs/{accession}.gz.log"
    params:
        extra="--skip-technical"
    threads: 6  # defaults to 6
    resources:
        mem_mb=4000,
        runtime=120,
    wrapper:
        "0.80.2/bio/sra-tools/fasterq-dump"


rule sample_sheet:
    input:
        "reads/{accession}_1.fastq.gz",
        "reads/{accession}_2.fastq.gz",
    output:
        "output/{accession}.csv"
    params:
        acc=lambda wildcards: wildcards.accession
    shell:
        """
        echo 'sample,group,short_reads_1,short_reads_2,long_reads' >> {output}
        echo '{params.acc},0,{input[0]},{input[1]},' >> {output}
        """


rule mag_pipeline:
    input:
        input="output/{accession}.csv",
    output:
        "results/Assembly/MEGAHIT/{accession}.contigs.fa.gz",
        "results/Assembly/MEGAHIT/{accession}.log",
    params:
        pipeline="nf-core/mag",
        revision="2.1.1",
        profile=["singularity"],
    handover: True
    wrapper:
        "https://raw.githubusercontent.com/hivlab/snakemake-wrappers/nf-profile/utils/nextflow"


rule trimmed_reads:
    input: 
        "results/Assembly/MEGAHIT/{accession}.log"
    output: 
        "results/trimmed_reads/{accession}.phix_removed.unmapped_1.fastq.gz",
        "results/trimmed_reads/{accession}.phix_removed.unmapped_2.fastq.gz",
    shell:
        """
        a=($(grep -o "(\/.*phix_removed.unmapped_[1,2].fastq.gz" {input[0]} |\
        sed "s/(//g" |\
        awk -F',' '{{ for(i=1;i<=NF;i++) print $i }}'))
        ln -sr "${{a[0]}}" {output[0]} \
        && ln -sr "${{a[1]}}" {output[1]}
        """


rule unzip:
    input: "results/Assembly/MEGAHIT/{accession}.contigs.fa.gz"
    output: "results/Assembly/MEGAHIT/{accession}.contigs.fa"
    shell:
        "zcat {input[0]} > {output[0]}"


rule metator:
    input:
        "results/trimmed_reads/{accession}.phix_removed.unmapped_1.fastq.gz",
        "results/trimmed_reads/{accession}.phix_removed.unmapped_2.fastq.gz",
        "results/Assembly/MEGAHIT/{accession}.contigs.fa",
    output:
        "results/metator/{accession}/bin_summary.txt",
        "results/metator/{accession}/contig_data_final.txt",
        "results/metator/{accession}/alignment_0.pairs",
    log:
        "logs/{accession}.metator.log"
    params:
        extra="--force",
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    container:
        "docker://koszullab/metator"
    threads: 8
    resources:
        mem_mb=44000,
        runtime=600,
    shell:
        """
        metator pipeline --forward='{input[0]}' --reverse='{input[1]}' --assembly='{input[2]}' --outdir='{params.outdir}' --threads={threads} {params.extra} 2> {log}
        """


rule virsorter2:
    input:
        "results/Assembly/MEGAHIT/{accession}.contigs.fa",
    output:
        "results/virsorter2/{accession}/final-viral-combined.fa",
        "results/virsorter2/{accession}/final-viral-score.tsv",
        "results/virsorter2/{accession}/final-viral-boundary.tsv",
    log:
        "logs/{accession}.virsorter2.log"
    params:
        extra="--min-score 0.5 --min-length 500",
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    container:
        "docker://jiarong/virsorter:latest"
    threads: 4
    resources:
        mem_mb=44000,
        runtime=600,
    shell:
        """
        virsorter run -w {params.outdir} {params.extra} -j {threads} all 2> {log}
        """


rule contactmap:
    input:
        "results/metator/{accession}/contig_data_final.txt",
        "results/metator/{accession}/alignment_0.pairs",
        "results/Assembly/MEGAHIT/{accession}.contigs.fa",
    output:
        directory("results/metator/{accession}/contact_map"),
    log:
        "logs/{accession}.metator.log"
    params:
        extra="",
        bin="MetaTOR_26_0",
        enzyme="MluCI",
    container:
        "docker://koszullab/metator"
    threads: 8
    resources:
        mem_mb=44000,
        runtime=600,
    shell:
        """
        metator contactmap -c {input[0]} -a {input[2]} -e {params.enzyme} -n {params.bin} -p {input[1]} -DfF -s 5000 -o {output[0]} --threads={threads} {params.extra} 2> {log}
        """

