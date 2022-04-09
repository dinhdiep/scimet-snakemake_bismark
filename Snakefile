'''
Snakefile
scimet-snakemake

Analysis of bisulfite reads from sciMETv2

2022-02-26 Copyright Dinh Diep

'''
import json
from os.path import join, basename, dirname
from subprocess import check_output
from itertools import chain

# Globals ---------------------------------------------------------------------

# Configurations for this run is indicated in this file.
configfile: 'config.yml'

# Full path to bismark reference folder.
BISMARK_REF = config['BISMARK_REF']

# Full path to folder with bowtie2 executables.
BOWTIE2_DIR = config['BOWTIE2_DIR']

min_map_score = config['MAP_SCORE']

# Samples and their corresponding filenames.
FILES = json.load(open(config['SAMPLES_JSON']))
SAMPLES = sorted(FILES.keys())


# Functions -------------------------------------------------------------------



# Rules -----------------------------------------------------------------------

rule all:
    input:
        'mapping_complete.txt'


checkpoint split:
    input:
    output:
        chunks=directory("processed_fastqs/{sample}")
    threads: 1
    run:
        # splitting files makes it possible to re-start jobs when it needs to be terminated in the middle of the step.
        r1 = FILES[wildcards.sample]['R1']
        r2 = FILES[wildcards.sample]['R2']
        shell("mkdir -p processed_fastqs/{wildcards.sample}")
        shell("mkdir -p status/{wildcards.sample}") 
        shell("mkdir -p bam_files/{wildcards.sample}") 
        shell('mkdir -p trimmed_fastqs/{wildcards.sample}')
        shell('zcat {r1} | split -l 400000000 - processed_fastqs/{wildcards.sample}/split.R1_')
        shell('zcat {r2} | split -l 400000000 - processed_fastqs/{wildcards.sample}/split.R2_')

 
rule trim:
    input:
        "processed_fastqs/{sample}/split.R1_{chunk}",
        "processed_fastqs/{sample}/split.R2_{chunk}"
    output:
        "trimmed_fastqs/{sample}/chunk_{chunk}.trimmed.R1.fq.gz",
        "trimmed_fastqs/{sample}/chunk_{chunk}.trimmed.R2.fq.gz"
    threads: 8
    run:
        shell('pigz -p {threads} {input[0]}')
        shell('pigz -p {threads} {input[1]}')
        shell('perl ./sciMETv2/sciMET_trim.pl -t {threads} -1 {input[0]}.gz -2 {input[1]}.gz -O trimmed_fastqs/{wildcards.sample}/chunk_{wildcards.chunk}')
        shell('rm {input[0]}.gz {input[1]}.gz')
        # make placeholders
        shell('echo "Done trimming" > {input[0]}')


rule align:
    input:
        "trimmed_fastqs/{sample}/chunk_{chunk}.trimmed.R1.fq.gz",
        "trimmed_fastqs/{sample}/chunk_{chunk}.trimmed.R2.fq.gz"
    output:
        "bam_files/{sample}/{sample}.chunk_{chunk}.nsrt.bam"
    threads: 8
    run:
       shell('perl ./sciMETv2/sciMET_align.pl -X -w {BOWTIE2_DIR} -t {threads} -R {BISMARK_REF} -1 {input[0]} -2 {input[1]} -O {wildcards.sample}.chunk_{wildcards.chunk}')
       shell('mv {wildcards.sample}.chunk_{wildcards.chunk}.nsrt.bam bam_files/{wildcards.sample}/')


def nsrt_bam_input(wildcards):
    checkpoint_output = checkpoints.split.get(**wildcards).output[0]
    return expand("bam_files/{sample}/{sample}.chunk_{chunk}.nsrt.bam",
           sample = wildcards.sample,
           chunk = glob_wildcards(os.path.join(checkpoint_output, "split.R1_{chunk}")).chunk)


rule merge_nsrt_bam:
    input:
        nsrt_bam_input
    output:
        "bam_files/{sample}/{sample}.merged.nsrt.bam",
    threads: 24
    run:
        shell('samtools merge -@ {threads} -n -c -p {output[0]} {input}')


rule rmdup:
    input:
        "bam_files/{sample}/{sample}.merged.nsrt.bam"
    output:
        "bam_files/{sample}/{sample}.mapping.complete"
    threads: 24
    run:
        shell('perl ./sciMETv2/sciMET_rmdup.pl -t {threads} -q {min_map_score} -O bam_files/{wildcards.sample}/{wildcards.sample} {input}')
        shell('echo {input} > {output}')


rule mapping_complete:
    input:
        e = expand("bam_files/{sample}/{sample}.mapping.complete", sample = SAMPLES)
    output:
        protected('mapping_complete.txt')
    run:
        shell('echo {input} > {output}')

