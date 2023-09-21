#!/bin/bash
# FASTQ FORMAT
#     @SEQ_ID
#     GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
#     +
#     !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65

fasta=$1
fastq="$(echo ${fasta} | sed 's/.fa$/.fq/g')"

LENGTH=150

seq_id=$(head -1 ${fasta} | sed 's/ /_/g')
seq=$(tail -n+2 ${fasta} | tr -d '\n' | head -c ${LENGTH})
qual=$(for ((i=1; i<=150; i++)); do echo -n "~"; done)

printf "${seq_id}\n${seq}\n+\n${qual}" > ${fastq}

echo "${fastq}"
