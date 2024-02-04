#!/bin/bash

data=/home/rccuser/shared/emily_files/S_nigra_norrna/rcorrector_output/
rsem_reference_output_file=/home/rccuser/shared/emily_files/RSEM_nigra/TrinityRefNigra.idx.fa

for fq in ${data}$1*_1.cor.fq.gz
        do
        base=$(basename $fq _1.cor.fq.gz)
        echo "Running RSEM for ${base} ..."
        time rsem-calculate-expression -p 16 --paired-end --bowtie2 \
                $fq ${data}${base}_2.cor.fq.gz \
                /home/rccuser/shared/emily_files/RSEM_nigra/TrinityRefNigra.idx.fa ${base}
done
