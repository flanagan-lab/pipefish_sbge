#!/bin/bash

#Create the arguments
ref_trans=$1 #Transcriptome .fasta file
index_dir=$2 #dir/basename the bt2 files will be written to
read_dir=$3 #location of the processed reads that will be aligned to index
sam_dir=$4 #Where you want the output SAM files to be stored
bam_dir=$5 #Where you want the output BAM files to be stored

#Build an index based on the transcriptome
echo "Indexing the refrence transcriptome"
time bowtie2-build --threads 16 $ref_trans $index_dir


#Map the reads back to the index
echo "Now beginning alignments ..."
for fq in $read_dir*_R1.fq.gz
	do

	#Extract sample name from the file
	sample=$(basename $fq _R1.fq.gz)

	#Echo the name of the sample that is currently running
	echo "Aligning reads for $sample"

	#Align this pair of reads
	time bowtie2 -x $index_dir \
		-1 $read_dir${sample}_R1.fq.gz \
		-2 $read_dir${sample}_R2.fq.gz \
		-S $sam_dir${sample}.sam --threads 16

	#Convert the SAM files to BAM files and coordinate sort
	echo "Converting sam file to bam file for ${sample}"
	time samtools view -T $ref_trans -b $sam_dir${sample}.sam > $bam_dir${sample}.bam
	echo "Sorting bam file for ${sample}"
	time samtools sort -@ 16 -o $bam_dir${sample}_sorted.bam $bam_dir${sample}.bam

done
