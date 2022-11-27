#!/bin/bash

trinity_output_file=$1
RSEM_output=$2

rsem-prepare-reference --bowtie2 $trinity_output_file $RSEM_output
