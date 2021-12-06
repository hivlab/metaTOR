
FASTQ=["SRR13435231"]

# Wrappers Github repo: https://github.com/avilab/virome-wrappers
WRAPPER_PREFIX = "https://raw.githubusercontent.com/avilab/virome-wrappers"

rule all:
    input: expand(["reads/{accession}_1.fastq.gz", "reads/{accession}_2.fastq.gz", "output/{accession}/contigs_fixed.fa", "output/{accession}/metator"], accession=FASTQ)


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

rule assembly:
    input:
        pe1 = "reads/{accession}_1.fastq.gz",
        pe2 = "reads/{accession}_2.fastq.gz",
    output:
        contigs="output/{accession}/final.contigs.fa",
    params:
        extra=(
            lambda wildcards, resources: f"--presets meta-large --min-contig-len 1000 --verbose -m {resources.mem_mb * 1048576}"
        ),
    threads: 8
    log:
        "logs/{accession}_assembly.log",
    shadow:
        "minimal"
    resources:
        runtime=lambda wildcards, input: round(2600 + 0.06 * input.size_mb),
        mem_mb=lambda wildcards, input: round(20000 + 2.22 * input.size_mb),
    wrapper:
        f"{WRAPPER_PREFIX}/v0.8.0/assembly/megahit"


rule fix_fasta:
    input:
        rules.assembly.output.contigs,
    output:
        "output/{accession}/contigs_fixed.fa",
    params:
        lambda wildcards: wildcards.accession,
    conda:
        f"{WRAPPER_PREFIX}/v0.8.0/subset_fasta/environment.yaml"
    resources:
        mem_mb=4000,
        runtime=120,
    script:
        "scripts/fix_fasta.py"


rule metator:
    input:
        "reads/{accession}_1.fastq.gz",
        "reads/{accession}_2.fastq.gz",
        "output/{accession}/final.contigs.fa",
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