#!/bin/bash
# Script to go from a TSV file created from ./staxis_to_ncbi_assemlby.sh containing 
# taxid, organism, and fasta download to the BLASTDB files needed to run blastn

# INPUT_FILE is the file w/ format from staxis_to_ncbi_assemlby.sh  - taxid, organism, refseq_type, assembly, url
# $ ./staxis_to_ncbi_assemlby.sh 1897061
# 1897061 Cryobacterium sp. SO1   na      Scaffold        https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/210/215/GCA_004210215.1_ASM421021v1/GCA_004210215.1_ASM421021v1_genomic.fna.gz
INPUT_FILE=$1

grep "https" ${INPUT_FILE} | \
  while IFS= read line; do
    staxid=$(echo "${line}" | cut -f1)
    species=$(echo "${line}" | cut -f2 | sed 's/ /_/g' | sed 's/\.//g')
    refseq_type=$(echo "${line}" | cut -f3 | sed 's/ /_/g')
    assembly=$(echo "${line}" | cut -f4 | sed 's/ /_/g')
    url=$(echo "${line}" | cut -f5 | sed 's/ /_/g')
    pushd . > /dev/null

    workdir="${staxid}___${species}___${refseq_type}___${assembly}"
    echo ${workdir}
    mkdir -p ${workdir}

    cd ${workdir}

    log_dir=$(realpath logs)
    mkdir -p ${log_dir}

    fasta_dir=$(realpath "fasta")
    mkdir -p ${fasta_dir}
    readme="README.md"
    echo "staxid=${staxid}" > ${readme}
    echo "species=${species}" >> ${readme}
    echo "refseq_type=${refseq_type}" >> ${readme}
    echo "assembly=${assembly}" >> ${readme}
    echo "" >> ${readme}
    echo "fasta download=${url}" >> ${readme}
    echo "" >> ${readme}

    fna_file=$(find . -type f -name "*.fna*")
    if [[ -f ${fna_file} ]]; then
      printf "\tdownloaded. Skipping download...\n"
    else
      cd ${fasta_dir}
      printf "\tdownloading ${url}\n"
      wget ${url} 2> ${log_dir}/log_download_${staxid}.out
      cd - > /dev/null
    fi

    input_blastdb_fasta_file="${fasta_dir}/${species}__${staxid}.fa"
    if [[ -f ${input_blastdb_fasta_file} ]]; then
      printf "\tFound: ${input_blastdb_fasta_file}. Skipping fasta...\n"
    else 
      printf "\tpreparing taxid fasta: $(basename ${input_blastdb_fasta_file})\n"
      gz_file=$(find . -type f -name "*.gz")
      gzcat ${gz_file} | \
        while IFS= read line; do
          if [[ ! -z $(echo "${line}" | grep -E "^>") ]]; then
            echo "${line}" | cut -d' ' -f1 >> ${input_blastdb_fasta_file}
          else
            echo "${line}" >> ${input_blastdb_fasta_file}
          fi
        done
      rm ${gz_file}
    fi

    blast_db_folder=$(realpath "blast_db")
    mkdir -p ${blast_db_folder}
    tax_id_map="${blast_db_folder}/taxid_map__${species}__${staxid}.txt"
    printf "\tPreparing tax ID map: $(basename ${tax_id_map})\n"
    grep ">" ${input_blastdb_fasta_file} | sed "s/^>//g" | sed "s/$/ ${staxid}/g" >> ${tax_id_map}

    blast_db_title="${species}__${staxid}"
    printf "\tCreating BLAST DB: ${blast_db_title}\n"
    CMD="makeblastdb -in ${input_blastdb_fasta_file} -parse_seqids -taxid_map ${tax_id_map} -title '${blast_db_title}' -dbtype nucl"
    printf "\t${CMD}\n"
    echo "${CMD}" >> ${readme}

    blast_db_log="${log_dir}/log_blastdb_creation___${blast_db_title}.out"
    eval ${CMD} > ${blast_db_log} 2>&1
    if [[ 0 -eq $? ]]; then
      printf "\tSUCCESS TAXID=${staxid} SPECIES=${species}\n"
    else
      printf "\tFAIL TAXID=${staxid} SPECIES=${species}\n\t\tLOGS=${blast_db_log}\n"
    fi

    echo "Moving BLASTDB files to separate directory"
    mv "${input_blastdb_fasta_file}.*" ${blast_db_folder}
    popd > /dev/null
  done

echo "Done." 
