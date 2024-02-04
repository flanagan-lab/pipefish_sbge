#!/bin/bash

#Create the arguements
input_dir=$1 #This should be the location of the RAW reads to be trimmed

for fq in $1*_R1_001.fastq.gz
        do
        base=$(basename $fq _R1_001.fastq.gz)
        echo "Running trimmomatic for ${base}..."
        time trimmomatic PE -threads 16 $fq $1${base}_R2_001.fastq.gz \
                $1trimmed/${base}_paired_R1_001.fastq.gz $1trimmed/${base}_unpaired_R1_001.fastq.gz \
                $1trimmed/${base}_paired_R2_001.fastq.gz $1trimmed/${base}_unpaired__R2_001.fastq.gz \
                HEADCROP:12 LEADING:3 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:50
done
