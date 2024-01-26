#!/bin/bash

#Create arguements
samples_file=$1 #File containing sample names and locations
out_dir_name=$2 #Desired name for the output directory

sudo docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq Trinity --seqType fq --max_memory 188G \
	--samples_file $1 \
	--CPU 16 \
	--output `pwd`/$2
