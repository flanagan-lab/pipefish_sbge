#!/bin/bash

#Create the arguements
input_dir_TR=$1 #Location of the .txt files that contain the trinity gene IDs
subset_fasta=$2 #Path to the subset_fasta_file script
assembly_file=$3 #Name/location of the de novo assembly
blast_database=$4
output_dir=$5 #Desired output directory for .fasta and BLAST files

## Loop through all the Trinity gene ID .txt files
for file in $input_dir_TR/*TRgenes.txt
	do

	#Extract the sample name from the file name
	sample=$(basename $file .txt)

	#Get the corresponding sequences for the Trinity Gene IDs
	echo "Running subset_fasta_file for ${sample}..."
	$subset_fasta -c -f $assembly_file -l $input_dir_TR/${sample}.txt

	#Rename the automated output from the subset_fasta_file script
	fasta_out=$(basename $assembly_file fasta)
	mv $fasta_out.subset.fasta ${sample}.fasta

	#Blast the sequences
	echo "Running BLAST for ${sample}..."
	blastn -db $blast_database -query ${sample}.fasta -out ${sample}_blast.txt \
		-evalue 0.001 \
		-num_threads 12 \
		-outfmt "6 qseqid qstart qend stitle sstart send evalue bitscore length pident gaps"

	#Move the outputs into desired output directory
	mv ${sample}* $output_dir

done
