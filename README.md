# short-read-seq_snakemake
SnakeMake pipeline for alignment, filtering, and peak-calling for short-read sequencing data.

## Description
This pipeline uses SnakeMake to automate processing of short-read sequencing datasets. The key steps and software packages are:

1. Quality control of sequencing reads (ChIPQC)
2. Sequencing adapter trimming (cta)
3. Alignment to reference genome (BWA mem)
4. Alignment quality evaluation (PICARD)
5. Alignment filtering (SAMTools)
6. Peak calling/signal generation (FSeq2)

The pipeline can be customized through the SnakeMake configuration file to apply to paired-end or single-end data, and to use any available genome assembly.

## Requirements
The pipeline requires the following software to be preinstalled:

Python >=2.7, and the following software packages:

1. fastqc
2. cta (<https://github.com/ParkerLab/cta>)
3. BWA (<http://bio-bwa.sourceforge.net/>)
4. picard (<https://broadinstitute.github.io/picard/>)
5. samtools (<http://www.htslib.org/>)
6. Fseq2 (<https://github.com/Boyle-Lab/F-Seq2>)
7. bedtools (<https://bedtools.readthedocs.io/en/latest/index.html>)


## Running the Pipeline

The SnakeMake pipeline relies on a configuration file in json format that specifies the source data and various run parameters, including the reference genome, blacklist (if any), etc. See the example below. Note that some sections are optional, and it is only necessary to provide one library, although multiple libraries may be included if desired.

```bash
{
    "blacklist": {  # (Optional) For each genome, a list of blacklisted regions in bed format
		    # (Not required by the pipeline, but should be used if they are available!!).
		    # These are used for peak filtering, and by ataqv.
        "hg19": [
            "/path/to/blacklist/wgEncodeDukeMapabilityRegionsExcludable.bed.gz", 
            "/path/to/blacklist/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
        ]
    },
    "bwa_index": { # (Required) path to BWA indices for each genome needed
        "hg19": "/path/to/bwa/index/hg19",
        "mm9": "/path/to/bwa/index/mm9"
    }, 
    "results": "/lab/work/porchard/atacseq", # (Optional) Path to the directory in which results should be placed (default is current working directory is used)
    "libraries": { # (Required) this is where the information for each library is given
        "100474___2156": { # unique ID for first library
            "genome": "hg19", # genome for first library
            "readgroups": { # readgroups for first library
			    # if the library was sequenced across several lanes, multiple readgroups
			    # can be provided, and they will be merged after mapping and before duplicate marking/filtering
			    # in this case, the library was sequenced across four lanes so four readgroups are provided.
                "100474___L1___2156": [ # list of the 2 fastq files for the first lane
                    "/path/to/data/fastq/100474_L001.1.fastq.gz", 
                    "/path/to/data/fastq/100474_L001.2.fastq.gz"
                ], 
                "100474___L2___2156": [
                    "/path/to/data/fastq/100474_L002.1.fastq.gz", 
                    "/path/to/data/fastq/100474_L002.2.fastq.gz"
                ]
            }
        }, 
        "100477___2156": { # second library begins here
            "genome": "mm9", 
            "readgroups": {
                "100477___L1___2156": [
                    "/path/to/data/fastq/100477_L001.1.fastq.gz", 
                    "/path/to/data/fastq/100477_L001.2.fastq.gz"
                ], 
                "100477___L2___2156": [
                    "/path/to/data/fastq/100477_L002.1.fastq.gz", 
                    "/path/to//data/fastq/100477_L002.2.fastq.gz"
                ]
            }
        }, 
    }
}

```

## Credits
This pipeline is based on the original Snakemake ATAC-seq pipeline developed by the Parker lab. The original version can be found here:
<https://github.com/ParkerLab/ATACseq-Snakemake>

