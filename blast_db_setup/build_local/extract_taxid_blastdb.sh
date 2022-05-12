#!/bin/bash
# Writes the BLAST DB of the input taxonomic ID
#   - Extracts only the sequences 
#   - Using those sequences, creates an isolated blastDB just for that taxonomic ID

help_string="\t./extract_taxid_blastdb.sh -t \${TAXID} [-s \${SPECIES}] [-d \${EXTRACT_DB}] [-h]\n"

while getopts ":t:s:d:h" opt; do
    case $opt in
        t) TAXID=${OPTARG}
        ;;
        s) SPECIES=${OPTARG}
        ;;
        d) EXTRACT_DB=${OPTARG}
        ;;
        h) printf "${help_string}" && exit 0  
        ;;
    esac 
done

if [[ -z ${TAXID} ]]; then
  echo "Please provide the taxonomic ID. Usage below,"
  printf "${help_string}"
  exit 1
fi
which blastdbcmd > /dev/null
if [[ $? -ne 0 ]]; then
  echo "PATH does not have blastdbcmd. Exiting"
  exit 1
fi
if [[ -z ${EXTRACT_DB} ]]; then
  EXTRACT_DB="nt"  # The nt BLAST DB is the default
fi
if [[ -z ${SPECIES} ]]; then
  SPECIES=$(blastdbcmd -taxids ${TAXID} -db ${EXTRACT_DB} -outfmt "%S" 2> /dev/null | head -1)
  if [[ $? -ne 0 ]]; then
    echo "[ERROR] Cannot determine species from TAXID=${TAXID}. Plesae verify taxonomic ID is correct, or provide the species"
    exit 1
  fi
fi

SPECIES="$(echo ${SPECIES} | sed 's/ /_/g' | sed 's/\.//g' | sed 's/\//_/g' | sed 's/(/_/g' | sed 's/)/_/g' | sed 's/\[/_/g' | sed 's/\]/_/g')"

echo "TAXID=${TAXID}"
echo "SPECIES=${SPECIES}"
echo "EXTRACT_DB=${EXTRACT_DB}"

pushd . > /dev/null

workdir="${TAXID}___${SPECIES}___${EXTRACT_DB}"
echo ${workdir}
mkdir -p ${workdir}

cd ${workdir}
fasta_dir=$(realpath fasta)
blast_db_folder=$(realpath "blast_db")

log_dir=$(realpath logs)
mkdir -p ${log_dir}

readme=$(realpath README.md)
echo "TAXID=${TAXID}" > ${readme}
echo "SPECIES=${SPECIES}" >> ${readme}
echo "EXTRACT_DB=${EXTRACT_DB}" >> ${readme}
echo "" >> ${readme}

mkdir -p ${fasta_dir}
input_blastdb_fasta_file="${fasta_dir}/${TAXID}___${SPECIES}.fa"
if [[ -f ${input_blastdb_fasta_file} ]]; then
  printf "\tFound: ${input_blastdb_fasta_file}. Skipping fasta...\n"
else 
  fasta_log="${log_dir}/fasta___${TAXID}___${SPECIES}.out"
  printf "\tpreparing taxid fasta: $(basename ${input_blastdb_fasta_file})\n"

  temp="${input_blastdb_fasta_file}.raw"
  FASTA_CMD="blastdbcmd -taxids ${TAXID} -db ${EXTRACT_DB} > ${temp}"
  echo "${FASTA_CMD}" >> ${readme}
  eval ${FASTA_CMD} > ${fasta_log} 2>&1
  if [[ 0 -ne $? ]]; then
    echo "[ERROR] FAILED fasta extraction from ${EXTRACT_DB}. See ${fasta_log}"
    rm -rf ${fasta_dir}
    exit 1
  fi

  printf "\tRemoving spaces from fasta description lines\n"
  FORMAT_CMD="sed 's/ .*//g' ${temp} > ${input_blastdb_fasta_file}"
  echo "${FASTA_CMD}" >> ${readme}
  eval ${FORMAT_CMD} >> ${fasta_log} 2>&1
  rm ${temp}
  if [[ 0 -ne $? ]]; then
    echo "[ERROR] FASTA formatting failed. See ${fasta_log}"
    rm -rf ${fasta_dir}
    exit 1
  fi
fi

tax_id_map="${fasta_dir}/taxid_map__${SPECIES}__${TAXID}.txt"
printf "\tPreparing tax ID map: $(basename ${tax_id_map})\n"
grep ">" ${input_blastdb_fasta_file} | sed "s/^>//g" | sed "s/$/\t${TAXID}/g" > ${tax_id_map}

mkdir -p ${blast_db_folder}
blast_db_title="${SPECIES}__${TAXID}"
cd ${blast_db_folder} 
printf "\tCreating BLAST DB: ${blast_db_title}\n"
CMD="makeblastdb -in ${input_blastdb_fasta_file} -parse_seqids -taxid_map ${tax_id_map} -title '${blast_db_title}' -dbtype nucl -out ${blast_db_title}"
printf "\t${CMD}\n"
echo "${CMD}" >> ${readme}
blast_db_log="${log_dir}/log_blastdb_creation___${blast_db_title}.out"
eval ${CMD} > ${blast_db_log} 2>&1
if [[ 0 -eq $? ]]; then
  printf "\tSUCCESS TAXID=${TAXID} SPECIES=${SPECIES}\n"
  printf "\tRemoving FASTA\n"
  rm ${input_blastdb_fasta_file}
else
  printf "\tFAIL TAXID=${TAXID} SPECIES=${SPECIES}\n\t\tLOGS=${blast_db_log}\n"
  rm -rf ${blast_db_folder}
fi
popd > /dev/null