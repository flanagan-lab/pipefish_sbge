#!/bin/bash

data=/home/rccuser/20220902_mRNASeq_PE150/trimmed_paired/

rsem-generate-data-matrix ${data}FUG10M2_paired.genes.results ${data}FUG11M2_paired.genes.results ${data}FUG11M4_paired.genes.results ${data}FUG12M1_paired.genes.results ${data}FUG15M5_paired.genes.results \
	${data}FUL10M2_paired.genes.results ${data}FUL11M2_paired.genes.results ${data}FUL11M4_paired.genes.results ${data}FUL12M1_paired.genes.results ${data}FUL15M5_paired.genes.results \
	${data}FUT10M2_paired.genes.results ${data}FUT11M2_paired.genes.results ${data}FUT11M4_paired.genes.results ${data}FUT12M1_paired.genes.results ${data}FUT15M5_paired.genes.results \
	${data}FUG11F1_paired.genes.results ${data}FUG13F1_paired.genes.results ${data}FUG13F4_paired.genes.results ${data}FUG2F2_paired.genes.results ${data}FUG3F2_paired.genes.results \
	${data}FUL11F1_paired.genes.results ${data}FUL13F1_paired.genes.results ${data}FUL13F4_paired.genes.results ${data}FUL2F2_paired.genes.results ${data}FUL3F2_paired.genes.results \
	${data}FUO11F1_paired.genes.results ${data}FUO13F1_paired.genes.results ${data}FUO13F4_paired.genes.results ${data}FUO2F2_paired.genes.results ${data}FUO3F2_paired.genes.results \
	> RSEM_fuscus_digi.counts.matrix
