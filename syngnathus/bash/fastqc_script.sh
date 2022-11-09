#!/bin/bash

data=/home/rccuser/20220902_mRNASeq_PE150/

for fq in ${data}trimmed_paired/*.gz
	do
	base=$(basename $fq)
	echo "Running fastqc for ${base} ..."
	time fastqc $fq -t 16 -o ${data}FastQC_trimmed
done
