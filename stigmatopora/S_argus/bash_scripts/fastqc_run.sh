#!/bin/bash

input_dir=$1 #location of the reads
output_dir=$2 #name of the desired output directory

for fq in $1/*.gz
        do
        base=$(basename $fq)
        echo "Running fastqc for ${base} ..."
        time fastqc $fq -t 16 -o $2
done
