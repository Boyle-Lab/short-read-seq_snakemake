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

## Dependencies
The pipeline requires the following software to be preinstalled:

Python >=2.7, and the following software packages:

1. fastqc
2. cta (<https://github.com/ParkerLab/cta>)
3. BWA (<http://bio-bwa.sourceforge.net/>)
4. picard (<https://broadinstitute.github.io/picard/>)
5. samtools (<http://www.htslib.org/>)
6. Fseq2 (<https://github.com/Boyle-Lab/F-Seq2>)
7. bedtools (<https://bedtools.readthedocs.io/en/latest/index.html>)

We strongly recommend installing the required packages in an Anaconda environment using Python >=3.7.


## Running the Pipeline

The SnakeMake pipeline relies on a configuration file in json format that specifies the source data and various run parameters, including the reference genome, blacklist (if any), etc. See the example below. Note that some sections are optional. Note that the pipeline expects fastq files for reads 1 and 2 to be designated with the _R1_001 and _R2_002 suffixes. This behavior can be modified by editing the wildcard values for the input and output specifications for the "trim" rule. 

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
    "results": "/path/to/results", # (Optional) Path to the directory in which results should be placed (default is current working directory is used)
    "libraries": { # (Required) this is where the information for the sequencing library is given
        "WT_Sperm_IP_with_MS-H3k9me2#2": { # unique ID for first library
            "genome": "hg19", # genome for first library
            "readgroups": { # readgroups for first library
			    # if the library was sequenced across several lanes, multiple readgroups
			    # can be provided, and they will be merged after mapping and before duplicate marking/filtering.
                "1619-SS-16": [
                    "/path/to/data/fastq/1619-SS-16_CGGACAAC-AATCCGGA_S20_R1_001.fastq.gz",
                    "/path/to/data/fastq/1619-SS-16_CGGACAAC-AATCCGGA_S20_R2_001.fastq.gz"
                ]
            }
        },
	"WT_Sperm_Input": {
            "genome": "mm10",
            "readgroups": {
                "1619-SS-21": [
                    "/path/to/data/fastq/1619-SS-21_ATGGCATG-GGTACCTT_S26_R1_001.fastq.gz",
                    "/path/to/data/fastq/1619-SS-21_ATGGCATG-GGTACCTT_S26_R2_001.fastq.gz"
                ]
            }
        }

    }
}
```

**IMPORTANT**: the basename for each fastq file must be unique.

To run the pipeline, you simply run snakemake with appropriate arguments indicating the name of the config file, location of the Snakefile, and number of cores. It is a good idea to wrap the command with nohup, and to redirect output to a log file, as shown in the example below. Note that the -p option to snakemake causes all shell commands to be printed as they are executed, which is useful for troubleshooting.

```bash
nohup snakemake -p --cores 8 --configfile snakemake_config.json --snakefile /path/to/Snakefile &> snakemake.nohup &
```

### build_snakefile.sh
To make the process of assembling library data into the JSON config file, we supply a shell script, build_snakefile.sh, which can be found in the scripts directory. This script is run as follows:

```bash
bash scripts/build_snakefile.sh sampleInfo.tsv /path/to/fastq/data /path/to/results/directory genome /path/to/bwa/index snakemake_config.json

```

You can substitute any filename you want for 'snakemake_config.json'

Example files are given given in `examples/`.

In case you are running on a cluster and need a cluster config file for
Snakemake, a template cluster config can be found in `src/` as well.

## Credits
This pipeline is based on the original Snakemake ATAC-seq pipeline developed by the Parker lab. The original version can be found here:
<https://github.com/ParkerLab/ATACseq-Snakemake>

