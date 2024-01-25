#!/bin/bash

#Create arguments
input_dir=$1
index_file=$2
output_dir=$3

for fq in $input_dir*_R1.fq.gz
	do

	#Extract sample name from the file
	sample=$(basename $fq _R1.fq.gz)

	#Echo sample name that is currently running
	echo "Processing sample ${sample} ..."

	#Quantify this pair of reads
	/usr/local/bin/salmon quant -i $index_file -l A \
		-1 $input_dir${sample}_R1.fq.gz \
		-2 $input_dir${sample}_R2.fq.gz \
		-p 16 -o $output_dir${sample}_quant

done

