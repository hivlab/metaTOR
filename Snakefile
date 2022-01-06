import os
from snakemake.utils import min_version
min_version("6.2")


report: "report/workflow.rst"

pepfile: "pep.yaml"
RUNS=pep.sample_table["run"].unique()
GROUPS=list(set(pep.sample_table.index.values))


localrules: all, sample_sheet, trimmed_reads, unzip
rule all:
    input: 
        expand(
            ["reads/{accession}_1.fastq.gz", 
            "reads/{accession}_2.fastq.gz"], 
            accession=RUNS
            ),
        expand(
            ["results/{group}.csv",
            "results/Assembly/MEGAHIT/{group}.contigs.fa.gz"],
            group=GROUPS
        )


rule get_fastq_pe_gz:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        "reads/{accession}_1.fastq.gz",
        "reads/{accession}_2.fastq.gz",
    log:
        "logs/{accession}.fasterq-dump.log"
    params:
        extra="--skip-technical"
    threads: 6  # defaults to 6
    resources:
        mem_mb=lambda wildcards, input: int(4000 + input.size / 1e6),
        runtime=120,
    wrapper:
        "0.80.2/bio/sra-tools/fasterq-dump"


rule sample_sheet:
    input:
        lambda wildcards: expand("reads/{accession}_{pair}.fastq.gz", accession=pep.sample_table.loc[[wildcards.group], "run"].tolist(), pair=[1,2])
    output:
        "results/{group}.csv"
    params:
        group=lambda wildcards: wildcards.group
    run:
        df=pep.sample_table.loc[[params.group], ["sample", "group", "short_reads_1", "short_reads_2", "long_reads"]]
        df.to_csv(output[0], index=False)


rule mag_pipeline:
    input:
        input="results/{group}.csv",
    output:
        "results/Assembly/MEGAHIT/{group}.contigs.fa.gz",
        "results/Assembly/MEGAHIT/{group}.log",
    params:
        pipeline="nf-core/mag",
        revision="2.1.1",
        profile=["singularity"],
        min_contig_size=1000,
        coassemble_group=True,
    handover: True
    wrapper:
        "https://raw.githubusercontent.com/hivlab/snakemake-wrappers/nf-profile/utils/nextflow"


rule trimmed_reads:
    input: 
        "results/Assembly/MEGAHIT/{group}.log"
    output: 
        "results/trimmed_reads/{run}.phix_removed.unmapped_1.fastq.gz",
        "results/trimmed_reads/{run}.phix_removed.unmapped_2.fastq.gz",
    shell:
        """
        a=($(grep -o "(\/.*phix_removed.unmapped_[1,2].fastq.gz" {input[0]} |\
        sed "s/(//g" |\
        awk -F',' '{{ for(i=1;i<=NF;i++) print $i }}'))
        ln -sr "${{a[0]}}" {output[0]} \
        && ln -sr "${{a[1]}}" {output[1]}
        """


rule unzip:
    input: "results/Assembly/MEGAHIT/{run}.contigs.fa.gz"
    output: "results/Assembly/MEGAHIT/{run}.contigs.fa"
    shell:
        "zcat {input[0]} > {output[0]}"


rule metator:
    input:
        "results/trimmed_reads/{run}.phix_removed.unmapped_1.fastq.gz",
        "results/trimmed_reads/{run}.phix_removed.unmapped_2.fastq.gz",
        "results/Assembly/MEGAHIT/{run}.contigs.fa",
    output:
        "results/metator/{run}/bin_summary.txt",
        "results/metator/{run}/contig_data_final.txt",
        "results/metator/{run}/alignment_0.pairs",
    log:
        "logs/{run}.metator.log"
    params:
        extra="--force",
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    container:
        "docker://koszullab/metator:latest"
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
        "results/Assembly/MEGAHIT/{run}.contigs.fa",
    output:
        "results/virsorter2/{run}/final-viral-combined.fa",
        "results/virsorter2/{run}/final-viral-score.tsv",
        "results/virsorter2/{run}/final-viral-boundary.tsv",
    log:
        "logs/{run}.virsorter2.log"
    params:
        extra="--min-score 0.5 --min-length 500",
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    container:
        "docker://jiarong/virsorter:2.2.3"
    threads: 4
    resources:
        mem_mb=44000,
        runtime=600,
    shell:
        """
        virsorter run -i {input[0]} -w {params.outdir} {params.extra} -j {threads} all 2> {log}
        """


rule contactmap:
    input:
        "results/metator/{run}/contig_data_final.txt",
        "results/metator/{run}/alignment_0.pairs",
        "results/Assembly/MEGAHIT/{run}.contigs.fa",
    output:
        directory("results/metator/{run}/contact_map"),
    log:
        "logs/{run}.metator.log"
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

