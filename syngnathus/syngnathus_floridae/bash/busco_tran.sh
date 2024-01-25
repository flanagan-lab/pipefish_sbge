#!/bin/bash

#Create arguments
transcriptome=$1 #Output fasta file from Trinity
lineage=$2 #chosen dataset for assessment
out_dir_name=$3 #Desired name for the output directory

busco -i $1 -l $2 -m tran -o $3 -c 16
