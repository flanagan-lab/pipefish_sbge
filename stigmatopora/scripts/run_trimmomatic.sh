#!/bin/bash

data=/mnt/BigData/OG7629/fastq ##This is the location of the raw reads
outloc=/mnt/BigData/OG7629/trimmed

for fq in ${data}/*_R1_001.fastq.gz
do
    echo $fq
    base=$(basename $fq _R1_001.fastq.gz)
    echo "Running trimmomatic for ${base}..."
    java -jar ~/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 6 $fq ${data}/${base}_R2_001.fastq.gz \
        ${outloc}/${base}_paired_R1.fastq.gz ${outloc}/${base}_unpaired_R1.fastq.gz \
        ${outloc}/${base}_paired_R2.fastq.gz ${outloc}/${base}_unpaired_R2.fastq.gz \
        ILLUMINACLIP:~/Programs/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 \
        HEADCROP:12 LEADING:3 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:50
done
