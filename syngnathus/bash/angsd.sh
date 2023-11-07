#!/bin/bash

#Choose what to run
ALLSFS=false
THETAS=false
INDSFS=true
MFFST=true

#Creating arguments
bam_dir=$1 #Folder where all of the desired BAM files are stored
angsd_dir=$2 #Path to the angsd executable
ref_trans=$3 #Reference transcriptome file (.fasta)
sfs_all_out=$4 #Desired name for the output SFS across all samples
stat_out=$5 #Desired name for the output of the Theta and/or Fst calculation

#Create the list of all SORTED bam files, just the female bam files and just the male bam files
ls ${bam_dir}*_sorted* > all_bams.txt
ls ${bam_dir}*M*_sorted* > male_bams.txt
ls ${bam_dir}FL*F*_sorted* > fem_bams.txt

if [ $ALLSFS = true ]; then
	#Estimate the allele frequencies from the BAM files
	echo "Estimating the site frequency spectrum across ALL samples"
	${angsd_dir}angsd -b all_bams.txt \
			-doSaf 1 \
			-anc ${ref_trans} \
			-minMapQ 20 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
			-minQ 13 -minInd 4 \
			-GL 1 -P 16 -out ${sfs_all_out}

	${angsd_dir}misc/realSFS ${sfs_all_out}.saf.idx -fold 1 -maxIter 100 -P 16 > ${sfs_all_out}.sfs

fi


if [ $THETAS = true ]; then
	#Calculating Thetas
	echo "Calculating thetas ..."
	${angsd_dir}misc/realSFS saf2theta ${sfs_all_out}.saf.idx \
			-sfs ${sfs_all_out}.sfs \
			-outname ${theta_out}

	#Calculating Tajima's D
	echo "Calculating Tajima's D..."
	${angsd_dir}misc/thetaStat do_stat ${stat_out}.thetas.idx

fi


if [ $INDSFS = true ]; then
	#Calculate the site frequency spectrum for males and females seperately
	echo "Estimating seperate SFS..."

	#First calculate per pop saf for each population (in our cases the different sexes)
	${angsd_dir}angsd -b fem_bams.txt -doSaf 1 -anc ${ref_trans} \
			-minMapQ 20 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
			-minQ 13 -minInd 7 \
			-GL 1 -P 12 -out ${sfs_all_out}_fem

	${angsd_dir}angsd -b male_bams.txt -doSaf 1 -anc ${ref_trans} \
			-minMapQ 20 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
			-minQ 13 -minInd 7 \
			-GL 1 -P 12 -out ${sfs_all_out}_mal

	#Calculate the 2dsfs prior
	${angsd_dir}misc/realSFS ${sfs_all_out}_fem.saf.idx ${sfs_all_out}_mal.saf.idx -fold 1 > fem.mal.ml

fi

if [ $MFFST = true ]; then
	#Calculate M-F Fsts
	echo "Calculating male-female Fsts..."

	#Prepare the Fst for easy window analysis etc.
	${angsd_dir}misc/realSFS fst index ${sfs_all_out}_fem.saf.idx ${sfs_all_out}_mal.saf.idx -sfs fem.mal.ml -fstout ${stat_out}_fm_fst

	#Get the global estimate
	${angsd_dir}misc/realSFS fst stats ${stat_out}_fm_fst.fst.idx

fi
