#!/bin/bash
ALIGNMENTS="/mnt/BigData/genomes/orthofinder/OrthoFinder/Results_Sep16/Single_Copy_Orthologue_Alignments"
OUTDIR="/mnt/BigData/genomes/orthofinder/OrthoFinder/Results_Sep16/Single_Copy_Orthologue_absrel2"
for aln in $ALIGNMENTS/*_aln.fa
do
  og=$(basename $aln _aln.fa)
  hyphy absrel --alignment $aln --tree SpeciesTree_rooted.txt  --code Universal --branches FG --output ${OUTDIR}/${og}.absrel
done
