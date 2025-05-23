#!/bin/bash

#Create arguments
ref_fasta=$1 #desired reference database
input_dir=$2 #Input directory with the location of the reads
output_dir=$3 #Desired location of the output

## Loop through all pairs of reads in the input directory
for pair in $2/*_fwd.fq.gz
        do

        #Extract the sample name from the file name
        sample=$(basename $pair _fwd.fq.gz)

        ##Echo the sample name it is currently running
        echo "Running Kraken2 for ${sample}..."

        #Define the paths to the input and output files for this sample
        input1=$2/${sample}_fwd.fq.gz
        input2=$2/${sample}_rev.fq.gz

        #Run Kraken2 on this pair of reads
        time kraken2 --threads 16 --db $1 --paired $input1 $input2 --unclassified-out $3/${sample}#.fq --report $3/${sample}.log

done
