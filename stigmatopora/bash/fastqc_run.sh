#!/bin/bash

data=/home/rccuser/shared/emily_files/S_nigra_raw/trimmed/paired/

for fq in ${data}*.gz
	do
	base=$(basename $fq)
	echo "Running fastqc for ${base} ..."
	time fastqc $fq -t 16 -o ${data}FastQC_trimmed
done
