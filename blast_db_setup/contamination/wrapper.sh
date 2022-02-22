#!/bin/bash

query=$1 

if [[ ! -f ${query} ]]; then
  echo "query=${query} is not a valid file"
  exit 1
fi

log=log_blastn.out
ulimit -n 8192
ulimit -f -n > ${log}

result=results.out


echo "${query}" >> ${log}
wc -l ${query} >> ${log}

start=$(date +"%s")
echo "start=${start}"

echo "start=${start}" >> ${log}

./ncbi-blast-2.12.0+/bin/blastn \
  -db ref_euk_rep_genomes \
  -num_threads 7 \
  -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames sblastnames sskingdoms" \
  -query ${query} \
  -out ${result} >> ${log} 2>&1

end=$(date +"%s")
echo "end=${end}"

let "time = end - start"
echo "time=${time}"

printf "Time=${time}\n\tstart=${start}\n\tend=${end}\n" >> ${log}
