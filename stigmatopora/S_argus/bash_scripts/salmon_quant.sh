#!/bin/bash

#Create arguments
input_dir=$1
index_file=$2
output_dir=$3

for fq in $1/*_1.cor.fq
        do

        #Extract sample name from the file
        sample=$(basename $fq _1.cor.fq)

        #Echo sample name that is currently running
        echo "Processing sample ${sample} ..."

        #Quantify this pair of reads
        /usr/local/bin/salmon quant -i $index_file -l A \
                -1 $1/${sample}_1.cor.fq \
                -2 $1/${sample}_2.cor.fq \
                -p 16 -o $3/${sample}_quant

done
