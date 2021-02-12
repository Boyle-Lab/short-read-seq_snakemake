import itertools
import os
import sys
import functools


# GENERIC DATA

ORGANISMS = {
    'rn4': 'rat',
    'rn5': 'rat',
    'mm9': 'mouse',
    'mm10': 'mouse',
    'hg19': 'human',
    'hg38': 'human'
}

AUTOSOMAL_REFERENCES = {
    'hg19': ['chr{}'.format(i) for i in range(1, 23)],
    'hg38': ['chr{}'.format(i) for i in range(1, 23)],
    'mm9': ['chr{}'.format(i) for i in range(1, 20)],
    'mm10': ['chr{}'.format(i) for i in range(1, 20)],
    'rn4': ['chr{}'.format(i) for i in range(1, 21)],
    'rn5': ['chr{}'.format(i) for i in range(1, 21)]
}

MACS2_GENOME_SIZE = {
    'rn4': 'mm',
    'rn5': 'mm',
    'mm9': 'mm',
    'mm10': 'mm',
    'hg19': 'hs',
    'hg38': 'hs'
}


# RESULT PATHS

prefix_results = functools.partial(os.path.join, config['results'])
FASTQC_DIR = prefix_results('fastqc')
TRIM_DIR = prefix_results('trim')
BWA_DIR = prefix_results('bwa')
MERGE_DIR = prefix_results('merge_readgroups')
MD_DIR = prefix_results('mark_duplicates')
PRUNE_DIR = prefix_results('prune')
MACS2_DIR = prefix_results('macs2')
LOG_DIR = prefix_results('logs')
VERSION_DIR = prefix_results('versions')


# Helper functions

def iterate_all_libraries():
    for library in sorted(config['libraries'].keys()):
        yield library


def iterate_library_readgroups(library):
    for rg in sorted(config['libraries'][library]['readgroups'].keys()):
        yield rg


def readgroup_to_library(readgroup):
    for library in iterate_all_libraries():
        for library_readgroup in iterate_library_readgroups(library):
            if readgroup == library_readgroup:
                return library


def iterate_all_readgroups():
    for library in iterate_all_libraries():
        for readgroup in iterate_library_readgroups(library):
            yield readgroup


def list_readgroup_fastqs(readgroup):
    library = readgroup_to_library(readgroup)
    return config['libraries'][library]['readgroups'][readgroup]


def iterate_all_fastqs():
    for readgroup in iterate_all_readgroups():
        for fastq in list_readgroup_fastqs(readgroup):
            yield fastq


def fastq_basename_to_fastq(fastq_basename):
    for fastq in iterate_all_fastqs():
        if fastq_basename == os.path.basename(fastq):
            return fastq
    sys.stderr.write('ERROR: could not find fastq corresponding to {}; exiting.\n'.format(fastq_basename))
    sys.exit()


def fastq_to_trimmed_fastq(fastq):
    trimmed_fastq_basename = os.path.basename(fastq).replace('.fastq.gz', '.trimmed.fastq.gz')
    return os.path.join(TRIM_DIR, trimmed_fastq_basename)


def get_genome(library):
    return config['libraries'][library]['genome']


def get_organism(genome):
    return ORGANISMS[genome]


def get_bwa_index(genome):
    return config['bwa_index'][genome]


def get_autosomes(genome):
    return AUTOSOMAL_REFERENCES[genome]

def get_all_chroms(genome):
    """ Get all chromosomes, excluding chrM. """
    tmp = AUTOSOMAL_REFERENCES[genome]
    tmp.extend(["chrX", "chrY"])
    return(tmp)


# Now the pipeline itself

rule all:
    input:
        # fastqc output
        [os.path.join(FASTQC_DIR, '{}_fastqc.zip'.format(os.path.basename(fastq).replace('.fastq.gz', ''))) for fastq in iterate_all_fastqs()],
	[os.path.join(PRUNE_DIR, '{}.pruned.bam.bai'.format(x)) for x in config['libraries']],
        # software versions
        os.path.join(VERSION_DIR, 'fastqc_version.txt'),
        os.path.join(VERSION_DIR, 'cta_version.txt'),
        os.path.join(VERSION_DIR, 'bwa_version.txt'),
        os.path.join(VERSION_DIR, 'macs2_version.txt'),
        os.path.join(VERSION_DIR, 'picard_version.txt'),
        os.path.join(VERSION_DIR, 'samtools_version.txt')


rule fastqc:
    input:
        lambda wildcards: fastq_basename_to_fastq('{}.fastq.gz'.format(wildcards.fastq_basename))
    output:
        os.path.join(FASTQC_DIR, '{fastq_basename}_fastqc.zip')
    params:
        outdir=FASTQC_DIR
    log:
        os.path.join(LOG_DIR, 'fastqc.{fastq_basename}.log')
    shell:
        'fastqc {input} -o {params.outdir} &> {log}'


rule trim:
    input:
        first = lambda wildcards: fastq_basename_to_fastq('{}_R1_001.fastq.gz'.format(wildcards.fastq_basename)),
        second = lambda wildcards: fastq_basename_to_fastq('{}_R2_001.fastq.gz'.format(wildcards.fastq_basename))
    output:
        first = os.path.join(TRIM_DIR, '{fastq_basename}_R1_001.trimmed.fastq.gz'),
        second = os.path.join(TRIM_DIR, '{fastq_basename}_R2_001.trimmed.fastq.gz')
    shell:
        'cta {input.first} {input.second} {output.first} {output.second}'


rule map:
    input:
        first = lambda wildcards: fastq_to_trimmed_fastq(list_readgroup_fastqs(wildcards.readgroup)[0]),
        second = lambda wildcards: fastq_to_trimmed_fastq(list_readgroup_fastqs(wildcards.readgroup)[1]),
        fasta = lambda wildcards: get_bwa_index(get_genome(wildcards.library))
    output:
        os.path.join(BWA_DIR, '{library}______{readgroup}.bam')
    params:
        sort_tmp = os.path.join(BWA_DIR, '{library}______{readgroup}.sort.tmp'),
        rg = '\\t'.join(['@RG', 'ID:{}'.format('{readgroup}'), 'LB:{}'.format('{library}')]),
	threads_bwa = workflow.cores * 0.5, threads_samtools = workflow.cores * 0.5
    threads: workflow.cores
    log:
        bwa = os.path.join(LOG_DIR, 'map.bwa.{library}______{readgroup}.log'),
        samtools = os.path.join(LOG_DIR, 'map.samtools.{library}______{readgroup}.log')
    shell:
        """bwa mem -M -R \'{params.rg}\' -t {params.threads_bwa} {input.fasta} {input.first} {input.second} 2> {log.bwa} | samtools sort -m 1g -@ {params.threads_samtools} -O bam -T {params.sort_tmp} -o {output} - 2> {log.samtools}"""


rule merge:
    input:
        lambda wildcards: [os.path.join(BWA_DIR, '{}______{}.bam'.format(wildcards.library, readgroup)) for readgroup in iterate_library_readgroups(wildcards.library)]
    output:
        bam = os.path.join(MERGE_DIR, '{library}.bam'),
        bam_index = os.path.join(MERGE_DIR, '{library}.bam.bai'),
    params:
        unsorted_merged_bam = os.path.join(MERGE_DIR, '{library}.unsorted.bam'),
        sort_tmp_prefix = os.path.join(MERGE_DIR, '{library}.sort')
    shell:
        """samtools merge {params.unsorted_merged_bam} {input}; ionice -c2 -n7 samtools sort -m 1G -O bam -T {params.sort_tmp_prefix} -o {output.bam} {params.unsorted_merged_bam}; samtools index {output.bam}; rm {params.unsorted_merged_bam}"""


rule mark_duplicates:
    input:
        bam = os.path.join(MERGE_DIR, '{library}.bam'),
        bam_index = os.path.join(MERGE_DIR, '{library}.bam.bai')
    output:
        bam = os.path.join(MD_DIR, '{library}.md.bam'),
        bam_index = os.path.join(MD_DIR, '{library}.md.bam.bai')
    params:
        metrics = os.path.join(MD_DIR, '{library}.metrics'),
        tmp_dir = MD_DIR
    shell:
        """java -jar $EBROOTPICARD/picard.jar MarkDuplicates I={input.bam} O={output.bam} ASSUME_SORTED=true METRICS_FILE={params.metrics} VALIDATION_STRINGENCY=LENIENT TMP_DIR={params.tmp_dir}; samtools index {output.bam}"""

rule prune:
    input:
        bam = os.path.join(MD_DIR, '{library}.md.bam'),
        bam_index = os.path.join(MD_DIR, '{library}.md.bam.bai')
    output:
        bam = os.path.join(PRUNE_DIR, '{library}.pruned.bam'),
        bam_index = os.path.join(PRUNE_DIR, '{library}.pruned.bam.bai')
    params:
        tmp_dir = MD_DIR,
        mapq = 30,
	chroms = lambda wildcards: get_all_chroms(get_genome('{}'.format(wildcards.library)))
    shell:
        """samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q {params.mapq} {input.bam} {params.autosomes} > {output.bam}; samtools index {output.bam}"""

rule linkbams:
    input:
        bam = os.path.join(PRUNE_DIR, '{library}.pruned.bam')
    output:
	bam = os.path.join(PRUNE_DIR, '{library}.pruned.bam')

rule versions:
    output:
        fastqc_version = os.path.join(VERSION_DIR, 'fastqc_version.txt'),
        cta_version = os.path.join(VERSION_DIR, 'cta_version.txt'),
        bwa_version = os.path.join(VERSION_DIR, 'bwa_version.txt'),
        picard_version = os.path.join(VERSION_DIR, 'picard_version.txt'),
        samtools_version = os.path.join(VERSION_DIR, 'samtools_version.txt'),
        macs2_version = os.path.join(VERSION_DIR, 'macs2_version.txt'),
    run:
        shell('fastqc --version &> {output.fastqc_version} || echo ""')
        shell('cta --version &> {output.cta_version} || echo ""')
        shell('bwa &> {output.bwa_version} || echo ""')
        shell('java -jar $EBROOTPICARD/picard.jar MarkDuplicates --version &> {output.picard_version} || echo ""')
        shell('samtools --version &> {output.samtools_version} || echo ""')
        shell('macs2 --version &> {output.macs2_version} || echo ""')


