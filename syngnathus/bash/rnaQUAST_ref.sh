#!/bin/bash

#Create arguements
rnaQUAST_path=$1 #Path to the rnaQuast.py file
txome_file=$2 #Trinity output .Trinity.fasta file
genome_ref=$3 #.fasta/.fa/.fna/.ffn/.frn file containg desired reference genome
gene_cord=$4 #file with gene coordinates in GTF/GFF format
out_dir=$5 #Name of the output directory

python $rnaQUAST_path --transcripts $txome_file \
	--reference $genome_ref \
	--gtf $gene_cord \
	--output_dir $out_dir --threads 16
