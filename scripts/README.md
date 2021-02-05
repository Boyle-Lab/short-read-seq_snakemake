# short-read-seq_snakemake supplementary scripts
Supplementary scripts for preparing config files and post-processing of data from the SnakeMake pipeline for alignment, filtering, and peak-calling for short-read sequencing data.

## Description
This folder contains helper scripts and a Makefile to automate peak-calling and downstream analysis of alignments from the main Snakemake pipeline. Steps performed include:

1. Peak calling with Fseq2.
2. Sample correlation matrix production with deeptools.
3. Local enrichment heatmap for +-2000bp relative to gene annotations with deeptools.


## Dependencies
The pipeline requires the following software to be preinstalled:

Python >=2.7, and the following software packages:

1. Fseq2 (<https://github.com/Boyle-Lab/F-Seq2>)
2. bedtools (<https://bedtools.readthedocs.io/en/latest/index.html>)
3. deeptools (</data/projects/protamine/data/mm10_refGene.bed>)

We strongly recommend installing the required packages in an Anaconda environment using Python >=3.7.


## Running the Pipeline

Directions for using build_snakefile.sh are given in the README.md for the main pipeline.

For all other steps, you need to edit the Makefile to point to the location of sampleInfo.tsv, pruned bam files produced by the main pipeline, and location of gene annotations file. These are in the first three lines of the Makefile.

To run all steps, run:
```bash
make
```

For peak calling without any further analysis, run:

```bash
make peaks
```

To create heatmaps and correlation matrices on existing	peaks, run:
```bash
make matrices
```

## Credits
All scripts written by Adam Diehl (adadiehl@umich.edu).