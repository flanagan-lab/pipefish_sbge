#!/bin/bash

# USAGE: Provide two arguments (in this order), 
# with paths relative to the project directory (above where the script is located): 
## 1. the name of the samples file (e.g., /nigra_gonads_trinity.txt)
## 2. the output directory (e.g., results/trinity_nigra_gonads/)

# Get the sample file
samples=$1
output=$2

# Setting the directory 
### move to the directory of the script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

cd ../
projdir=`pwd`


# run Trinity from within the docker
sudo docker run --rm -v/home/rccuser/:/home/rccuser/ trinityrnaseq/trinityrnaseq Trinity \
    --seqType fq \
    --normalize_by_read_set \
    --max_memory 62G \
    --samples_file ${projdir}/${samples} \
    --CPU 12 \
    --output ${projdir}/${output}

