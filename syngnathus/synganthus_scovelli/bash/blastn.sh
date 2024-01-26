#!/bin/bash

#Create arguements
blast_database=$1 #Database that is used as a reference to BLAST against
sequences_fasta_file=$2 #Sequences you want to BLAST
output_file_name=$3


blastn -db $blast_database \
	-query $sequences_fasta_file \
	-out $output_file_name \
	-evalue 0.001 \
	-num_threads 12 \
	-outfmt "6 qseqid qstart qend stitle sstart send evalue bitscore length pident gaps"
