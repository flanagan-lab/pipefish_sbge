### Creating a for loop to run multiple samples at a time
#!/bin/bash

#Create arguments
rcorrector_path=$1 #Path to the run_rcorrector.pl file
input_dir=$2 #Path to the input directory containing the reads
output_dir=$3 #Path to the desired output location

##Loop through all pairs of reads in the directory
for pair in $2/*_1.fq
        do

        #Extract sample name from the file name
        sample=$(basename $pair _1.fq)

        #Echo the sample name that is currently running
        echo "Running rcorrector for ${sample} ..."

        #Run rcorrector on this pair of reads
        time perl $1/run_rcorrector.pl -t 16 -1 $2/${sample}_1.fq -2 $2/${sample}_2.fq -od $3

done
