#!/bin/bash


INPUT="taxid_name.tsv"

tail -n+2 ${INPUT} | \
  while IFS= read line; do
    tax_id=$(echo "${line}" | cut -f1)
    name=$(echo "${line}" | cut -f2)
    family=$(echo "${line}" | cut -f3)
    fasta=$(./get_fastas.sh ${tax_id} "${name}")
    fasta_finder=$(ls ${tax_id}__*)
    if [[ ! -f ${fasta_finder} ]]; then
      continue
    else
      fastq=$(./make_fastqs.sh ${fasta_finder})
      printf "${tax_id}\t${name}\t${family}\t${fasta}\t${fastq}\n"
    fi
    

    
  done