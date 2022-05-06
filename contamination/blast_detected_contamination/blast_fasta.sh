#!/bin/bash

help_string="./blast_fasta.sh -f FASTA [-d DB] [-b blast_bin]"
options="\t-d, blast database: e.g. \"ref_euk_rep_genomes\", \"nt\"\n\t-b, blast binary: path to blastn binary\n"

while getopts ":f:d:b:h" opt; do
    case $opt in
        f) fasta_file=${OPTARG}
        ;;
        d) DB=${OPTARG}
        ;;
        b) BLASTN=${OPTARG}
        ;;
        h) printf "${help_string}" && exit 0
        ;;
    esac 
done

if [[ ! -f ${fasta_file} ]]; then
  echo "fasta_file=${fasta_file} is not a valid file"
  echo "${help_string}"
  printf "${options}"
  exit 1
fi
if [[ -z ${DB} ]]; then
  # Default to nt
  DB=nt
fi
if [[ -z ${BLASTN} ]]; then
  BLASTN="blastn"
fi

which ${BLASTN}
if [[ 0 -ne $? ]]; then
  echo "${BLASTN} is not valid. Please point to a valid blastn binary w/ \"-b /path/to/blastn\""
  exit 1
fi

echo "Blasting to DB=\"${DB}\""

log=log_blastn.out
ulimit -n 8192
ulimit -f -n > ${log}

result=blast_results.tsv

echo "${fasta_file}" >> ${log}
wc -l ${fasta_file} >> ${log}

start=$(date +"%s")
echo "start=${start}"

echo "start=${start}" >> ${log}

${BLASTN} \
  -db "${DB}" \
  -num_threads 7 \
  -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames sblastnames sskingdoms" \
  -query ${fasta_file} \
  -out ${result} >> ${log} 2>&1
if [[ 0 -ne $? ]]; then
  echo ""
  echo "[FAILED] See ${log}"
  echo ""
fi

end=$(date +"%s")
echo "end=${end}"

let "time = end - start"
echo "time=${time}"

printf "Time=${time}\n\tstart=${start}\n\tend=${end}\n" >> ${log}
