#!/bin/bash

trinity_gene_map=$1
salmon_quant_files=$2
out_file_name=$3

sudo docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq /bin/sh -c "cd /home/rccuser/shared/emily_files/ && /usr/local/bin/util/abundance_estimates_t>                --est_method salmon \
                --gene_trans_map $trinity_gene_map \
                --quant_files $salmon_quant_files \
                --out_prefix $out_file_name"
