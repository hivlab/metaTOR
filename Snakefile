from snakemake.utils import min_version
min_version("6.2")


FASTQ=["SRR13435231"]


localrules: all, sample_sheet
rule all:
    input: expand(["reads/{accession}_1.fastq.gz", "reads/{accession}_2.fastq.gz", "results/Assembly/MEGAHIT/{accession}.contigs.fa.gz", "output/{accession}/metator"], accession=FASTQ)


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


rule metator:
    input:
        "results/trimmed_reads/{accession}.phix_removed.unmapped_1.fastq.gz",
        "results/trimmed_reads/{accession}.phix_removed.unmapped_2.fastq.gz",
        "results/Assembly/MEGAHIT/{accession}.contigs.fa.gz",
    output:
        directory("output/{accession}/metator"),
    container:
        "docker://koszullab/metator"
    threads: 8
    resources:
        mem_mb=4000,
        runtime=120,
    shell:
        """
        metator pipeline -1 {input[0]} -2 {input[1]} -a {input[2]} -o {output} --threads {threads}
        """