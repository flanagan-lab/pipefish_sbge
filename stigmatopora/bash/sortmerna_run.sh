#!/bin/bash

input_dir=/home/rccuser/shared/emily_files/S_nigra_raw/trimmed/paired/
output_dir=/home/rccuser/shared/emily_files/S_nigra_raw/trimmed/paired/sortmerna/

ref_fasta=/home/rccuser/shared/emily_files/rRNA_databases_v4.3.4/smr_v4.3_fast_db.fasta

for pair in $input_dir/*_R1_001.fastq.gz
	do

	sample=$(basename $pair _R1_001.fastq.gz)

	echo "Running sortmerna for ${sample}..."

	rm -rf /home/rccuser/sortmerna

	input1=$input_dir/${sample}_R1_001.fastq.gz
	input2=$input_dir/${sample}_R2_001.fastq.gz

	
	output=$output_dir/${sample}/
		
	# mkdir $output_dir/${sample}/

	time sortmerna --ref $ref_fasta	--reads $input1 --reads $input2 --fastx --other --aligned --workdir $output --threads 16 --out2 

done
