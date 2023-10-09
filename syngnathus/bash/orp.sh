#!/bin/bash

#Create arguements
oyster_mk_file=$1 #FULL path to the oyster.mk file
read_1=$2
read_2=$3
output_name=$4

#Run oyster.mk
$oyster_mk_file TPM_FILT=1 \
	STRAND=RF \
	MEM=188 \
	CPU=16 \
	READ1=$read_1 \
	READ2=$read_2 \
	RUNOUT=$output_name \
	LINEAGE=actinopterygii_odb10
