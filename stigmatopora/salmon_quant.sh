#!/bin/bash

#Create arguments
input_dir=$1
index_file=$2
output_dir=$3

for fq in $1/*_1.cor.fq.gz
        do

        #Extract sample name from the file
        sample=$(basename $fq _1.cor.fq.gz)

        #Echo sample name that is currently running
        echo "Processing sample ${sample} ..."

        #Quantify this pair of reads
        salmon quant -i $2 -l A \
                -1 $1/${sample}_1.cor.fq.gz \
                -2 $1/${sample}_2.cor.fq.gz \
                -p 16 --validateMappings -o $3/${sample}_quant

done
