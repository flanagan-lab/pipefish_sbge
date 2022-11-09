#!/bin/bash

data=/home/rccuser/20220902_mRNASeq_PE150/

sudo docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq Trinity --seqType fq --max_memory 188G \
	--samples_file ${data}trimmed_paired/FU_trinity_samples.txt \
	--CPU 12 \
	--output `pwd`/trinity_out_dir_fuscus
