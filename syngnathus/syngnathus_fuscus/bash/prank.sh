#!/bin/bash

# paths relative to location of script
FASTA_DIR="/mnt/BigData/genomes/orthofinder/OrthoFinder/Results_Aug15/Single_Copy_Orthologue_Seqs_NT"
SPP_TREE="../OrthoFinder/Results_Aug15/Species_Tree/SpeciesTree_rooted.txt"
OUT_DIR="/mnt/BigData/genomes/orthofinder/OrthoFinder/Results_Aug15/Single_Copy_Orthologue_Alignments"

# move to location of script
cd "${0%/*}" 

# check if the output directory does not exist
if [ ! -d "$OUT_DIR" ]; then
  mkdir $OUT_DIR # make the dir if it doesn't already exist
fi

for fa in $FASTA_DIR/*fasta
do
	ID=$(basename $fa .fasta) # remove the extension
	prank -d=$fa -o=$OUT_DIR/$ID -translate -shortnames
done

