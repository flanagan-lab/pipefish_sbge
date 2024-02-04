#!/bin/bash

trinity_output_file=$1 ##the .Trinity.fasta file
RSEM_output=$2 ##Ex. TrinityRefNigra

rsem-prepare-reference --bowtie2 $trinity_output_file $RSEM_output
