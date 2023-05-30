#!/bin/bash

fnames="ERR031940_1.filt.fastq.gz ERR031940_2.filt.fastq.gz"

echo "RETRIEVING FQs"
for fname in ${fnames}; do
  if [[ ! -f ${fname} ]]; then
    echo "getting ${fname}"
    wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00284/sequence_read/${fname}  2> /dev/null
    echo "  done."
  else
    echo "${fname} already downloaded. Skipping"
  fi
done

