#!/bin/bash

#Create arguments
transcriptome_file=$1
desired_index_name=$2
kmer_length=$3

salmon index -t $1 -i $2 -k $3
