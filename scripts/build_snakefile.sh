#!/bin/bash

# Build the snakefile for alignment, filtering, and peak-calling using the short-read snakemake pipeline.

SAMPLES=$1
FASTQ_PATH=$2
RESULTS=$3
GENOME=$4
IDX_PATH=$5
OUT=$6

IFS=$'\n'

echo | awk -F '\t' -v genome=$GENOME -v idx_path=$IDX_PATH -v res_path=$RESULTS -v fastq_path=$FASTQ_PATH '{printf "{\n\t\"bwa_index\": {\n\t\t\"%s\": \"%s\"\n\t},\n\t\"results\": \"%s\",\n\t\"libraries\": {\n", genome, idx_path, res_path}' > $OUT

i=1
LINES=$(cat $SAMPLES | grep -v -w "^#" | wc -l)
for REC in $(cat $SAMPLES | grep -v -w "^#"); do
    read -d "\n" RG LIB SNAME <<< $(echo $REC | awk -F '\t' '{printf "%s\n%s\n%s", $1, $2, $3}')

    R1=$(ls $FASTQ_PATH/$RG\_*R1*)  # Finds the fastq file for read1
    R2=$(ls $FASTQ_PATH/$RG\_*R2*)  # Finds the fastq file for read2
    
    echo $REC | awk -F '\t' -v genome=$GENOME -v lines=$LINES -v CR=$i -v fq1=$R1 -v fq2=$R2 '{printf "\t\t\"%s\": {\n\t\t\t\"genome\": \"%s\",\n\t\t\t\"readgroups\": {\n\t\t\t\t\"%s\": [\n\t\t\t\t\t\"%s\",\n\t\t\t\t\t\"%s\"\n\t\t\t\t]\n\t\t\t}\n\t\t}", $2, genome, $1, fq1, fq2; if (CR < lines) {printf ",\n"} else {printf "\n"}}' >> $OUT
    i=$((i + 1))
done

printf "\t}\n}\n" >> $OUT
