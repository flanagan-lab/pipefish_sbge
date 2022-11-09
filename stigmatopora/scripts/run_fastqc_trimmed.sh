#!/bin/bash

dir=/mnt/BigData/OG7629/trimmed
qcdir=/mnt/BigData/OG7629/trimmed_qc

for fq in ${dir}/*paired*.gz #The location of the PAIRED trimmed reads
        do
        base=$(basename $fq)
        echo "Running fastqc for ${base} ..."
        ~/Programs/FastQC/fastqc $fq -t 2 -o ${qcdir}/ #Where I want the output files to be stored
done