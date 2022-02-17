#!/bin/bash

input_sam=$1

dc_reads=dc_reads.sam
while IFS= read line; do
  # skip all headers (beginning w/ @)
  if [[ -z $(echo ${line} | grep -E "^@") ]]; then
    ds_cigar=$(echo ${line} | cut -d' ' -f6 | grep -oE "^[0-9]+S.*S$")
    if [[ ! -z ${ds_cigar} ]]; then
      echo ${line} >> ${dc_reads}
    fi
  fi
done < ${input_sam}

