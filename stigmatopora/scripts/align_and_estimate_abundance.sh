#!/bin/bash

## Needs to run in sudo mode ##

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"
cd ../results/

## Determines which RSEM step to run!
STEP=$1

if [[ "$STEP" == 1 ]]; then
####### Prepare the reference ########
        docker run -v/home/rccuser/:/home/rccuser trinityrnaseq/trinityrnaseq /usr/local/bin/util/align_and_estimate_abundance.pl --transcripts trinity_nigra_gonads.Trinity.fasta --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference
fi


#
#   ## Run the alignment and abundance estimation (assumes reference has already been prepped, errors-out if prepped reference not located.)
#
#    /usr/local/bin/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left reads_1.fq --right reads_2.fq --est_method RSEM --aln_method bowtie --trinity_mode --output_dir rsem_outdir
#
##  ## prep the reference and run the alignment/estimation
#
#    /usr/local/bin/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left reads_1.fq --right reads_2.fq --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir rsem_outdir
#
#   ## Use a samples.txt file:
#
#    /usr/local/bin/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --est_method RSEM --aln_method bowtie2 --prep_reference --trinity_mode --samples_file samples.txt --seqType fq  
#
