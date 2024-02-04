#!/bin/bash

data=/home/rccuser/shared/emily_files/
trinity_fasta_file=S_nigra_trinity.Trinity.fasta
trinity_stats_output=trinity_stats_output

sudo docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq /usr/local/bin/util/TrinityStats.pl
        ${data}$trinity_fasta_file
        >> $trinity_stats_output
