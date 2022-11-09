#!/bin/bash

data=/home/rccuser/20220902_mRNASeq_PE150/

for fq in ${data}*_R1.fastq.gz
	do
	base=$(basename $fq _R1.fastq.gz)
	echo "Running trimmomatic for ${base}..."
	time trimmomatic PE -threads 16 $fq ${data}${base}_R2.fastq.gz \
		${data}trimmed/${base}_paired_R1.fastq.gz ${data}trimmed/${base}_unpaired_R1.fastq.gz \
		${data}trimmed/${base}_paired_R2.fastq.gz ${data}trimmed/${base}_unpaired_R2.fastq.gz \
		ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 HEADCROP:12 LEADING:3 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:50
done
