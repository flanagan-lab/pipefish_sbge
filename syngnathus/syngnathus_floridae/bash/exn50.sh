#!/bin/bash

isoform_expression_matrix=$1
trinity_fasta_output=$2
ExN50_output=$3

sudo docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
        $isoform_expression_matrix \
        $trinity_fasta_output | tee $ExN50_output
