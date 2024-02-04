#!/bin/bash

#Create arguments
input_dir=$1 #location of the reads
ref_fasta=$2 #desired reference fasta
output_dir_rrna=$3 #desired location for the reads that are rRNA
output_dir_norrna=$4 #desired location for the reads that are NOT rRNA

## Loop through all pairs of reads in the input directory
for pair in $1/*_R1_001.fastq.gz
        do

        #Extract the sample name from the file name
        sample=$(basename $pair _R1_001.fastq.gz)

        #Extract Fish ID
        ID=$(basename $pair _paired_R1_001.fastq.gz)

        ##Echo the sample name it is currently running
        echo "Running SortMeRNA for ${ID}..."

        #Define the paths to the input and output files for this sample
        input1=$1/${sample}_R1_001.fastq.gz
        input2=$1/${sample}_R2_001.fastq.gz

        #Run SortMeRNA on this pair of reads
        time sortmerna --threads 16 --ref $2 --reads $input1 --reads $input2 --fastx --aligned $3/${ID} --other $4/${ID} --out2

        #Remove SortMeRNA intermediate files before running again
        rm -r /home/rccuser/sortmerna/run/kvdb

done
