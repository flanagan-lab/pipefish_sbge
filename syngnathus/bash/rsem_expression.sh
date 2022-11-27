#!/bin/bash

data=/home/rccuser/20220902_mRNASeq_PE150/trimmed_paired/
FU-or_FL=$1
rsem_reference_output_file=$2

for fq in ${data}$1*_R1.fastq.gz
	do
	base=$(basename $fq _R1.fastq.gz)
	echo "Running RSEM for ${base}..."
	time rsem-calculate-expression -p 10 --paired-end --bowtie2 \
		$fq ${data}${base}_R2.fastq.gz \
		$2 ${base}
done
