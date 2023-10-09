#!/bin/bash

#Create arguements
rnaQUAST_path=$1 #Path to the rnaQuast.py file
txome_file=$2 #Trinity output .Trinity.fasta file
input_dir=$3
out_dir=$4 #Name of the output directory

for fq in $input_dir*_R1.fq.gz
	do

	#Extract the sample name from the file name
	sample=$(basename $fq _R1.fq.gz)

	##Echo the sample its currently running
	echo "Running rnaQUAST for ${sample}"

	#Define the paths to the input and output reads for this pair
	read1=$input_dir${sample}_R1.fq.gz
	read2=$input_dir${sample}_R2.fq.gz

	#Run rnaQUAST on the pair of reads
	python $rnaQUAST_path --transcripts $txome_file \
		--left_reads $read1 \
		--right_reads $read2 \
		--output_dir $out_dir${sample} --threads 16 --gene_mark

done
