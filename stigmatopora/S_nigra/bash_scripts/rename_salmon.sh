#!/bin/bash

#Create arguements
input_dir=$1 #This is the location of all the Salmon quant folders that were generated
out_dir1=$2 #Desired location for the new renamed quant.sf files
out_dir2=$3 #Desired location for the new renamed salmon_quant.log files

for dir in $input_dir*_quant
    do

    echo "$dir"

    #Extract the sample name from the directory
    sample=$(basename $dir _quant)

    #Echo the sample that is being processed
    echo "Rename sample ${sample} ..."

    #Rename the quant.sf file
    mv $dir/quant.sf $out_dir1/${sample}_quant.sf

    #Rename the log file
    mv $dir/logs/salmon_quant.log $out_dir2/${sample}_salmon_quant.log

done
