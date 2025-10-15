#!/bin/bash

ALIGNMENTS="/mnt/BigData/genomes/orthofinder/OrthoFinder/Results_Sep16/Single_Copy_Orthologue_Alignments"
ref=single_copy_DEG.csv
species=("H_comes" "H_zosterae" "S_acus" "S_scovelli" "S_typhle" "S_fuscus")

for fa in $ALIGNMENTS/*best.nuc.fas
do
  og=$(basename $fa .best.nuc.fas)
  
  # remove stop codons
  hyphy ~/Programs/hyphy-2.5.78/res/TemplateBatchFiles/CleanStopCodons.bf Universal $fa "Disallow stops" $ALIGNMENTS/${og}_aln.fa

  # get the gene names for each species 
  genes=()
  genes+=( $(grep $og $ref | perl -pe 's/^.*\,OG\d+\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*),(.*)/\1/g') )
  genes+=( $(grep $og $ref | perl -pe 's/^.*\,OG\d+\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*),(.*)/\2/g') )
  genes+=( $(grep $og $ref | perl -pe 's/^.*\,OG\d+\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*),(.*)/\3/g') )
  genes+=( $(grep $og $ref | perl -pe 's/^.*\,OG\d+\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*),(.*)/\4/g') )
  genes+=( $(grep $og $ref | perl -pe 's/^.*\,OG\d+\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*),(.*)/\5/g') )
  genes+=( $(grep $og $ref | perl -pe 's/^.*\,OG\d+\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*)\,(.*),(.*)/\6/g') )


  # create a file with the regex to change the file names
  echo "#!/bin/bash" > fix_names_tmp.sh
  # for each species name, match it to ghe correct gene
  for i in "${!genes[@]}"; do
  	spp=$( echo "${species[$i]}" )
	g=$( echo  "${genes[$i]}" )
	echo "sed -ie 's/^>.*$g.*$/>$spp/g' $ALIGNMENTS/${og}_aln.fa" >> fix_names_tmp.sh
  done
  # run the file to change the names
  ./fix_names_tmp.sh

done




	
