#!/bin/bash

#Activate the conda environment
source activate orp

#List the reads
R1=$(ls $HOME/docker/coley_files/fuscus_kmer_corrected/*R1*.fq.gz | perl -pe 's/\n/ /g' | perl -pe 's/\s$//g')
R2=$(ls $HOME/docker/coley_files/fuscus_kmer_corrected/*R2*.fq.gz | perl -pe 's/\n/ /g' | perl -pe 's/\s$//g')

#Create arguements
oyster_mk_file=$1 #FULL path to the oyster.mk file
output_name=$2

#Run oyster.mk
$oyster_mk_file TPM_FILT=1 \
	STRAND=RF \
	MEM=188 \
	CPU=16 \
	READ1=$R1 \
	READ2=$R2 \
	RUNOUT=$output_name \
	LINEAGE=actinopterygii_odb10



