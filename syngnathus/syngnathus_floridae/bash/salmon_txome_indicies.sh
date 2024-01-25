#!/bin/bash

#Create arguments
transcriptome_file=$1 #This is the output .Trinity.fasta file from the txome assembly
desired_index_name=$2
kmer_length=$3 #Recommended 31 on the Salmon user guide

/usr/local/bin/salmon index -t $transcriptome_file \
			-i $desired_index_name \
			-k $kmer_length -p 16
