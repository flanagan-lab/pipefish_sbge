#!/bin/bash

data=/home/rccuser/20220902_mRNASeq_PE150/trimmed_paired/

rsem-generate-data-matrix ${data}FLG3M5_paired.genes.results ${data}FLG3M7_paired.genes.results ${data}FLG3M8_paired.genes.results ${data}FLG4M3_paired.genes.results ${data}FLG4M4_paired.genes.results \
        ${data}FLL3M5_paired.genes.results ${data}FLL3M7_paired.genes.results ${data}FLL3M8_paired.genes.results ${data}FLL4M3_paired.genes.results ${data}FLL4M4_paired.genes.results \
        ${data}FLT2M3_paired.genes.results ${data}FLT3M5_paired.genes.results ${data}FLT4M4_paired.genes.results ${data}FLT5M3_paired.genes.results ${data}FLT8M7_paired.genes.results \
        ${data}FLG2F7_paired.genes.results ${data}FLG3F1_paired.genes.results ${data}FLG3F2_paired.genes.results ${data}FLG8F3_paired.genes.results \
        ${data}FLL2F7_paired.genes.results ${data}FLL3F1_paired.genes.results ${data}FLL3F2_paired.genes.results ${data}FLL3F4_paired.genes.results ${data}FLL8F3_paired.genes.results \
        ${data}FLO2F7_paired.genes.results ${data}FLO3F1_paired.genes.results ${data}FLO3F2_paired.genes.results ${data}FLO3F4_paired.genes.results ${data}FLO8F3_paired.genes.results \
        > RSEM_floridae_digi.counts.matrix
