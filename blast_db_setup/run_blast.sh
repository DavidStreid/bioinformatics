#!/bin/bash

help_string="./run_blast sample.fa [/path/to/BLASTDB]"

fasta=$1
blastdb_path=$2	# PATH TO BLASTDB to use
blastdb_name=$3

echo "INPUTS"
printf "\tQuery=${fasta}\n"
printf "\tBLASTDB=${blastdb_path}\n"
printf "\tDatabase=${blastdb_name}\n"

if [[ ! -f ${fasta} ]]; then
  echo "Please provide a fasta file"
  printf "\t${help_string}\n"
  exit 1
fi
base=$(basename ${fasta} | cut -d'.' -f1)
result="blast_results__${base}.tsv"
log="log_blast__${base}.out"
if [[ ! -z ${blastdb_path} ]]; then
  echo "exporting BLASTDB=${blastdb_path}"
  export BLASTDB=${blastdb_path}
fi
echo "BLASTDB=${BLASTDB}"
if [[ -d ${BLASTDB} ]]; then
  echo "BLASTDB=${BLASTDB}" > ${log}
else
  >&2 echo "[ERROR] BLASTDB is not set. Please see configuration parameters - https://www.ncbi.nlm.nih.gov/books/NBK569858/"
  printf "\t${help_string}\n"
  exit 1
fi
if [[ -z ${blastdb_name} ]]; then
  echo "No blastdb name provided (check \"ls \$BLASTDB\" for what is available). Using nt"
  blastdb_name=nt
fi



ulimit -n 8192
ulimit -f -n >> ${log}
echo "${fasta}"  >> ${log}

wc -l ${fasta} >> ${log}



echo "RUNNING BLAST..."

FORMAT="6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames sblastnames sskingdoms"

CMD="blastn -db ${blastdb_name} -num_threads 7 -query ${fasta} -outfmt \"${FORMAT}\" -out ${result}"
echo "${CMD}"
echo "logging to ${log}"
eval "${CMD}" >> ${log} 2>&1
if [[ $? -eq 0 ]]; then
  echo "SUCCESS"
else
  echo "FAIL - see ${log}"
fi

echo "Done."
