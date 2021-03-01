#!/bin/bash

# Creates symbolic links to fastq files following naming conventions used in the processing
# pipeline. Sample data are in SAMPLES in the root folder.

SAMPLES=$1
BAM_PATH=$2  # Use the fully-qulaified absolute path here!

IFS=$'\n'

# Run a loop over samples to create the links.
for REC in $(cat $SAMPLES | grep -v '^#' | tr -d '\r'); do
    read -d "\n" ID NAME <<< $(echo $REC | awk -F '\t' '{printf "%s\n%s", $3, $2}')

    R1=$(ls $BAM_PATH/$NAME.*.bam)
    #printf "%s\t%s\n" $ID $R1
    ln -s $R1 $BAM_PATH/$ID.pruned.bam

done
