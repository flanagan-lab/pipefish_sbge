#!/bin/bash

R1=$(ls *R1*.fq.gz | perl -pe 's/\n/,/g' | perl -pe 's/,$//g')
R2=$(ls *R2*.fq.gz | perl -pe 's/\n/,/g' | perl -pe 's/,$//g')

#Create arguments
transcriptome=$1 #Output fasta file from Trinity

transrate --assembly $1 --left $R1 --right $R2 --threads 16




