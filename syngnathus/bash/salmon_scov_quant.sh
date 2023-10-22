#!/bin/bash

#Create arguments
input_dir=$1
index_file=$2
output_dir=$3

for fq in $input_dir*_R1.fq.gz
do
	#If a file has the lane informtion do one thing and if not run as normal
	if [[ $fq == *"L00"* ]]; then
		#If the sample is Lane 1 (L001) then continue, if it is Lane 2 (L002) stop
		if [[ $fq == *"L001"* ]]; then
			#Extract sample name from the file
			sample=$(basename $fq _L001_R1.fq.gz)

			#Count the number of files with that sample name
			num_files=$(ls $input_dir${sample}*.fq.gz | wc -l)

			#If the number of files is greater than 2 list L001 and L002 samples for -1 and -2
			if [ $num_files -gt 2 ]; then
				#Echo sample name that is currently running
				echo "Processing sample ${sample} ..."

				#Quantify this pair of reads (including the files from both lanes)
				/usr/local/bin/salmon quant -i $index_file -l A \
					-1 $input_dir${sample}_L001_R1.fq.gz $input_dir${sample}_L002_R1.fq.gz \
					-2 $input_dir${sample}_L001_R2.fq.gz $input_dir${sample}_L002_R2.fq.gz \
					-p 16 -o $output_dir${sample}_quant
			else
				#Echo sample name that is currently running
				echo "Processing sample ${sample} ..."

				#Quantify this pair of reads (Only across one lane)
				/usr/local/bin/salmon quant -i $index_file -l A \
					-1 $input_dir${sample}_L001_R1.fq.gz \
					-2 $input_dir${sample}_L001_R2.fq.gz \
					-p 16 -o $output_dir${sample}_quant
			fi
		fi
	else
		#Run the samples as normal, extract and echo sample name from the file
		sample=$(basename $fq _R1.fq.gz)
		echo "Processing sample ${sample} ..."

		#Quantify this pair of reads
		/usr/local/bin/salmon quant -i $index_file -l A \
			-1 $input_dir${sample}_R1.fq.gz  \
			-2 $input_dir${sample}_R2.fq.gz \
			-p 16 -o $output_dir${sample}_quant
	fi
done

