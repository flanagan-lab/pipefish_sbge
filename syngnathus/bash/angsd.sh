#!/bin/bash

#Creating arguments
bam_dir=$1 #Folder where all of the desired BAM files are stored
angsd_dir=$2 #Path to the angsd executable
ref_trans=$3 #Reference transcriptome file (.fasta)
sfs_all_out=$4 #Desired name for the output SFS across all samples
theta_out=$5 #Desired name for the output of the Theta calculation

#Create the list of all SORTED bam files, just the female bam files and just the male bam files
ls ${bam_dir}*_sorted* > all_bams.txt
ls ${bam_dir}*M*_sorted* > male_bams.txt
ls ${bam_dir}FL*F*_sorted* > fem_bams.txt

#Estimate the allele frequencies from the BAM files
echo "Estimating the site frequency spectrum across ALL samples"
${angsd_dir}angsd -b all_bams.txt \
		-doSaf 1 \
		-anc ${ref_trans} \
		-minMapQ 20 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
		-minQ 13 -minInd 4 \
		-GL 1 -P 16 -out ${sfs_all_out}

${angsd_dir}misc/realSFS ${sfs_all_out}.saf.idx -fold 1 -maxIter 100 -P 16 > ${sfs_all_out}.sfs

#Calculating Thetas
echo "Calculating thetas ..."
${angsd_dir}misc/realSFS saf2theta ${sfs_all_out}.saf.idx \
		-sfs ${sfs_all_out}.sfs \
		-outname ${theta_out}

#Calculating Tajima's D
echo "Calculating Tajima's D..."
${angsd_dir}misc/thetaStat do_stat ${theta_out}.thetas.idx

