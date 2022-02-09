#!/bin/bash

available_versions=$(curl ftp://ftp.ncbi.nih.gov/blast/executables/blast+/ 2>/dev/null | rev | cut -d' ' -f1 | rev | grep -oE "[0-9]+.[0-9]+.[0-9]+")
latest_version=$(curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/ 2> /dev/null | grep LATEST | sed 's/.*-> //g')

db_name_options="\t+-----------------------------+------------------------------------------------+\n"
db_name_options+="\t File Name                    | Content Description                           \n"
db_name_options+="\t+-----------------------------+------------------------------------------------+\n"
db_name_options+="\tnr.*tar.gz                    | Non-redundant protein sequences from GenPept,\n"
db_name_options+="\t                                Swissprot, PIR, PDF, PDB, and NCBI RefSeq\n"
db_name_options+="\tnt.*tar.gz                    | Partially non-redundant nucleotide sequences from \n"
db_name_options+="\t                                all traditional divisions of GenBank, EMBL, and DDBJ\n"
db_name_options+="\t                                excluding GSS,STS, PAT, EST, HTG, and WGS.\n"
db_name_options+="\tlandmark.tar.gz               | Proteome of 27 model organisms, see \n"
db_name_options+="\t                                https://blast.ncbi.nlm.nih.gov/smartblast/smartBlast.cgi?CMD=Web&PAGE_TYPE=BlastDocs#searchSets\n"
db_name_options+="\t16S_ribosomal_RNA             | 16S ribosomal RNA (Bacteria and Archaea type strains)\n"
db_name_options+="\t18S_fungal_sequences.tar.gz   | 18S ribosomal RNA sequences (SSU) from Fungi type and reference material (BioProject PRJNA39195)\n"
db_name_options+="\t28S_fungal_sequences.tar.gz   | 28S ribosomal RNA sequences (LSU) from Fungi type and reference material (BioProject PRJNA51803)\n"
db_name_options+="\tITS_RefSeq_Fungi.tar.gz       | Internal transcribed spacer region (ITS) from Fungi type and reference material (BioProject PRJNA177353)\n"
db_name_options+="\tITS_eukaryote_sequences.tar.gz| Internal transcribed spacer region (ITS) for eukaryotic sequences\n"
db_name_options+="\tLSU_eukaryote_rRNA.tar.gz     | Large subunit ribosomal RNA sequences for eukaryotic sequences\n"
db_name_options+="\tLSU_prokaryote_rRNA.tar.gz    | Large subunit ribosomal RNA sequences for prokaryotic sequences\n"
db_name_options+="\tSSU_eukaryote_rRNA.tar.gz     | Small subunit ribosomal RNA sequences for eukaryotic sequences\n"
db_name_options+="\tref_euk_rep_genomes*tar.gz    | Refseq Representative Eukaryotic genomes (1000+ organisms)\n"
db_name_options+="\tref_prok_rep_genomes*tar.gz   | Refseq Representative Prokaryotic genomes (5700+ organisms)\n"
db_name_options+="\tref_viruses_rep_genomes*tar.gz   | Refseq Representative Virus genomes (9000+ organisms)\n"
db_name_options+="\tref_viroids_rep_genomes*tar.gz   | Refseq Representative Viroid genomes (46 organisms)\n"
db_name_options+="\trefseq_protein.*tar.gz        | NCBI protein reference sequences\n"
db_name_options+="\trefseq_rna.*tar.gz            | NCBI Transcript reference sequences\n"
db_name_options+="\tswissprot.tar.gz              | Swiss-Prot sequence database (last major update)\n"
db_name_options+="\tpataa.*tar.gz                 | Patent protein sequences\n"
db_name_options+="\tpatnt.*tar.gz                 | Patent nucleotide sequences. Both patent databases\n"
db_name_options+="\t                                are directly from the USPTO, or from the EPO/JPO\n"
db_name_options+="\t                                via EMBL/DDBJ\n"
db_name_options+="\tpdbaa.*tar.gz                 | Sequences for the protein structure from the\n"
db_name_options+="\t                                Protein Data Bank\n"
db_name_options+="\tpdbnt.*tar.gz                 | Sequences for the nucleotide structure from the \n"
db_name_options+="\t                                Protein Data Bank. They are NOT the protein coding\n"
db_name_options+="\t                                sequences for the corresponding pdbaa entries.\n"
db_name_options+="\t+-----------------------------+------------------------------------------------+\n"
db_name_options+="\tFor more info, see https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html"

version_string="\tAvailable blast+ Versions: $(echo ${available_versions})\n\tLatest: ${latest_version}\n"
help_string="./setup_blast.sh (-v <version>) (-o <os>) (-d <database_name>) (-o <output_path>)\n"
help_string+="\tversions: $(echo ${available_versions})\n\t\tlatest (default): ${latest_version}\n"
help_string+="\tos: win64, x64-linux, x64-macosx, x64-win64\n\t\tDefault: x64-linux\n"
help_string+="\toutput_path: directory to write database files\n"
help_string+="\tdatabase_names: See https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html\n"


while getopts ":v:o:d:p:h" opt; do
    case $opt in
        v) version=${OPTARG}
        ;;
        o) os_input=${OPTARG}
        ;;
        d) db_name=${OPTARG}
        ;;
        p) out_path=${OPTARG}
        ;;
        h) printf "${help_string}\nOr, try './setup_blast.sh' to get the latest executables and take it from there\n"; exit 0
        ;;
    esac 
done

if [[ -z ${version} || ${version} == "l" ]]; then
  download_version=${latest_version}
else
  download_version=${version}
fi

if [[ -z ${os_input} ]]; then
  os="x64-linux"
else
  os=${os_input}
fi

if [[ -z ${out_path} ]]; then
  out_path="./"
fi

has_version=$(echo "${available_versions}" | grep -oE "^${download_version}$")

if [[ -z ${has_version} ]]; then
  echo "[ERROR] Invalid version: ${download_version}"
  printf "${version_string}"
  exit 1
fi 
untarred_file="ncbi-blast-${download_version}+"
tar_file="${untarred_file}-${os}.tar.gz"
echo "blast+ Version: ${download_version}"
echo "TAR file: ${tar_file}"
echo "version=${version} (latest: ${latest_version})"
echo "os_input=${os}"
if [[ ! -z ${db_name} ]]; then
  if [[ -z $(echo "${db_name_options}" | grep "${db_name}") ]]; then
    echo "[ERROR] Invalid db_name. See below,"
    printf "${db_name_options}"
    echo ""
    exit 1
  fi
  echo "db_name=${db_name}"
  echo "out_path=${out_path}"
fi
echo ""

if [[ -f ${tar_file} ]]; then
  echo "Detected downloaded file... skipping download"
elif [[ -f ${untarred_file} ]]; then
  echo "Detected extracted file... skipping download"
else
  CMD="curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${download_version}/${tar_file} -o ${tar_file}"
  echo "Downloading..."
  printf "\t${CMD}\n"
  eval ${CMD} 2> /dev/null
  if [[ $? -ne 0 ]]; then
    printf "\nDownload failed. Check that version (-v) and os (-o) are valid, or don't specify them to use defaults\n\tReceived version=${download_version} os=${os}\n"
    printf "\tNote - older version don't have all version, e.g. x64-macosx is not present in 2.2.* versions\n"
    echo "Exiting..."
    exit 1
  else
    echo "Download complete."
  fi
fi

 
if [[ -d ${untarred_file} ]]; then
  echo "Nothing to untar - Using ${untarred_file}"
elif [[ -f ${tar_file} ]]; then
  echo "Extracting contents..."
  tar -zxvf ${tar_file} 2> /dev/null
  echo "Extraction done."
else
  echo "No tarred/untarred file. Exiting..."
  exit 1
fi

download_script=$(realpath ncbi-blast-2.12.0+/bin/update_blastdb.pl)

if [[ ! -z ${db_name} ]]; then
  echo "Downloading & extracting files for database='${db_name}'"
  echo "Writing to ${out_path}"
  CMD="perl ${download_script} ${db_name}"
  echo ${CMD}
  cd ${out_path}
  log_dir=logs
  mkdir ${log_dir}
  eval ${CMD} > ${log_dir}/download_${db_name}.out 2>&1
  files=$(find . -type f -name "*.tar.gz" | sort)
  total=$(echo ${files} | cut ' ' '\n' | wc -l | grep -oE [0-9])
  for f in ${files}; do
    index=$(echo ${f} | xargs basename | cut -d'.' -f2)
    echo "[${index}] Extracting ${f}"
    tar -zxvf ${f} -C . >> ${log_dir}/extract.out 2>&1
    if [[ $? -ne 0 ]]; then
      echo "Extracted ${index}/${total}"
      echo "Failed to extract ${f}"
      echo "Extract ${f} manually and all remaining '*.tar.gz' files"
      echo "[WARN] DO NOT RE-RUN THIS SCRIPT"
      exit 1
    fi
    printf "\t[SUCCESS] Removing ${f}\n"
    rm ${f}
  done
  cd -
  echo "Successful download & extraction"
  echo "Export directory of db files as BLASTDB and run blast executables"
  echo "Blast executables: $(realpath ncbi-blast-2.12.0+/bin)"
  echo "EXPORT BLASTDB=$(realpath ${out_dir})"
else
  echo ""
  echo "For more info, see https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html"
  echo ""
  echo "Use this script to download databases"
  echo "perl ${download_script} <DB_NAME>"
  echo "  e.g."
  echo "    perl ${download_script} ref_euk_rep_genomes 	# Downloads all eukaryotic genomes"
fi
