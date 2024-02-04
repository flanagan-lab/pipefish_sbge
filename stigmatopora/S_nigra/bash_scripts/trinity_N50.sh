#!/bin/bash

##Create the arguments
trinity_fasta_file=$1 #This is the output .fasta file from your Trinity run
trinity_stats_output=$2 #Desired name for your output results

sudo docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq /usr/local/bin/util/TrinityStats.pl \
        $trinity_fasta_file \
        >> $trinity_stats_output
