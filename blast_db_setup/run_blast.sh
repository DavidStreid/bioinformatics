#!/bin/bash

help_string="./run_blast sample.fa [/path/to/BLASTDB]"

fasta=$1
blastdb=$2	# PATH TO BLASTDB to use

if [[ ! -f ${fasta} ]]; then
  echo "Please provide a fasta file"
  printf "\t${help_string}\n"
  exit 1
fi

base=$(basename ${fasta} | cut -d'.' -f1)
result="blast_results__${base}.tsv"
log="log_blast__${base}.out"
if [[ ! -z ${blastdb} ]]; then
  echo "exporting BLASTDB=${blastdb}"
  export BLASTDB=${blastdb}
fi
if [[ -d ${BLASTDB} ]]; then
  echo "BLASTDB=${BLASTDB}" > ${log}
else
  >&2 echo "[ERROR] BLASTDB is not set. Please see configuration parameters - https://www.ncbi.nlm.nih.gov/books/NBK569858/"
  printf "\t${help_string}\n"
  exit 1
fi

echo "BLASTDB=${BLASTDB}"

ulimit -n 8192
ulimit -f -n >> ${log}
echo "${fasta}"  >> ${log}

wc -l ${fasta} >> ${log}

echo "INPUTS"
printf "\tBLASTDB=${BLASTDB}\n"
printf "\tQUERY=${fasta}\n"

echo "RUNNING BLAST..."
blastn \
  -db nt \
  -num_threads 7 \
  -query ${fasta} \
  -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames sblastnames sskingdoms" \
  -out ${result} >> ${log} 2>&1

echo "Done."
