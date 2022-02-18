#!/bin/bash

input_sam=$1

if [[ ! -f ${input_sam} ]]; then
  echo "Please pass a SAM file"
  exit 1
fi

dc_reads='double_clipped_reads.sam'
output_fasta='double_clipped_reads.fa'

echo "Creating ${dc_reads}"
while IFS= read line; do
  # skip all headers (beginning w/ @)
  if [[ -z $(echo ${line} | grep -E "^@") ]]; then
    ds_cigar=$(echo ${line} | cut -d' ' -f6 | grep -oE "^[0-9]+S.*S$")
    if [[ ! -z ${ds_cigar} ]]; then
      echo ${line} >> ${dc_reads}
    fi
  fi
done < ${input_sam}

echo "Creating ${output_fasta}"
while IFS= read line; do
  CIGAR=$(echo ${line} | cut -d' ' -f6)
  QNAME=$(echo ${line} | cut -d' ' -f1)
  ID=">${QNAME}___${CIGAR}"
  echo ${ID} >> ${output_fasta}
  echo ${line} | cut -d' ' -f10 >> ${output_fasta}
done < ${dc_reads}
