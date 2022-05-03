#!/bin/bash

fasta_file=$1 

if [[ ! -f ${fasta_file} ]]; then
  echo "fasta_file=${fasta_file} is not a valid file"
  exit 1
fi

DB=ref_euk_rep_genomes
echo "Blasting to DB=${DB}"

log=log_blastn.out
ulimit -n 8192
ulimit -f -n > ${log}

result=blast_results.tsv


echo "${fasta_file}" >> ${log}
wc -l ${fasta_file} >> ${log}

start=$(date +"%s")
echo "start=${start}"

echo "start=${start}" >> ${log}

./ncbi-blast-2.12.0+/bin/blastn \
  -db ${DB} \
  -num_threads 7 \
  -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames sblastnames sskingdoms" \
  -query ${fasta_file} \
  -out ${result} >> ${log} 2>&1

end=$(date +"%s")
echo "end=${end}"

let "time = end - start"
echo "time=${time}"

printf "Time=${time}\n\tstart=${start}\n\tend=${end}\n" >> ${log}
