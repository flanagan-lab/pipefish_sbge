#!/bin/bash

# USAGE: Provide two arguments (in this order), 
# with paths relative to location of this script: 
## 1. the name of the samples file (e.g., ../nigra_gonads_trinity.txt)
## 2. the output directory (e.g., ../results/trinity_nigra_gonads/)

# Get the sample file
samples=$1
output=$2

# Setting the directory 
### move to the directory of the script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR


# run Trinity from within the docker
echo "sudo docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq Trinity \
    --seqType fq \
    --max_memory 30G \
    --samples_file ${samples} \
    --CPU 12 \
    --output ${output}"