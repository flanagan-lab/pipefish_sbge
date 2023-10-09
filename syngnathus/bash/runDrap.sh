#!/bin/bash

#Generate list of reads(comma seperated)
R1=$(ls *R1*.fq.gz | perl -pe 's/\n/,\/home\/rccuser\/shared\/coley_files\/fuscus_kmer_corrected\//g' | perl -pe 's/,\/home\/rccuser\/shared\/coley_files\/fuscus_kmer_corrected\/$//g')
R2=$(ls *R2*.fq.gz | perl -pe 's/\n/,\/home\/rccuser\/shared\/coley_files\/fuscus_kmer_corrected\//g' | perl -pe 's/,\/home\/rccuser\/shared\/coley_files\/fuscus_kmer_corrected\/$//g')

#Generate arguements
out_dir_name=$1

sudo docker run --rm -v`pwd`:`pwd` sigenae/drap runDrap --no-trim --outdir `pwd`/$1 \
		--R1 `pwd`/$R1 \
		--R2 `pwd`/$R2 \
		--dbg trinity --norm-mem 188 --dbg-mem 188
