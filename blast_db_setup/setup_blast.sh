#!/bin/bash
#!/bin/bash

available_versions=$(curl ftp://ftp.ncbi.nih.gov/blast/executables/blast+/ 2>/dev/null | rev | cut -d' ' -f1 | rev | grep -oE "[0-9]+.[0-9]+.[0-9]+")
latest_version=$(curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/ 2> /dev/null | grep LATEST | sed 's/.*-> //g')

version_string="\tAvailable blast+ Versions: $(echo ${available_versions})\n\tLatest: ${latest_version}\n"
help_string="./setup_blast.sh (-v <version>) (-o <os>) (-d <database_name>) (-p <output_path>)\n"
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
        d) list_of_dbs_to_download=${OPTARG}
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
if [[ ! -z ${list_of_dbs_to_download} ]]; then
  available_databases=$(curl ftp://ftp.ncbi.nlm.nih.gov/blast/db/  2>/dev/null | rev | cut -d' ' -f1 | rev | cut -d'.' -f1 | sort | uniq)
  valid_db_names=""
  for db_name in ${list_of_dbs_to_download}; do
    if [[ -z $(echo "${available_databases}" | grep -E "^${db_name}$") ]]; then
      echo "[WARN] Invalid db_name: ${db_name}"
    else
      valid_db_names+="${db_name} "
    fi
  done

  if [[ -z ${valid_db_names} ]]; then
    echo "[WARN] Will not download any DBs"
    printf "\tInvalid DB parameter=\"${list_of_dbs_to_download}\"\n"
    echo "Valid Databases Below"
    echo "${available_databases}"
    echo ""
  else
    list_of_dbs_to_download=${valid_db_names}
    echo "QUEUED_DATABASE_DOWNLOADS=${list_of_dbs_to_download}"
  fi
  out_path=$(realpath ${out_path})
  echo "DB_DIR=${out_path}"
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

binary_dir=$(realpath ncbi-blast-${latest_version}+/bin)
download_script="${binary_dir}/update_blastdb.pl"

if [[ ! -z ${list_of_dbs_to_download} ]]; then
  log_dir=logs
  for db_name in ${list_of_dbs_to_download}; do
    echo "Downloading & extracting files for database='${db_name}'"
    echo "Writing to ${out_path}"
    CMD="perl ${download_script} ${db_name}"
    echo ${CMD}
    cd ${out_path}  
    mkdir -p ${log_dir}
    eval ${CMD} > ${log_dir}/download_${db_name}.out 2>&1
    files=$(find . -type f -name "${db_name}*.tar.gz" | sort)
    total=$(echo ${files} | sed 's/ /\n/g' | wc -l | grep -oE [0-9])
    for f in ${files}; do
      index=$(echo ${f} | xargs basename | cut -d'.' -f2)
      echo "[${index}] Extracting ${f}"
      tar -zxvf ${f} -C . >> ${log_dir}/extract.out 2>&1
      if [[ $? -ne 0 ]]; then
        echo "Extracted ${index}/${total}"
        echo "Failed to extract ${f}"
        echo "Download/Extract ${f} manually and all remaining '*.tar.gz' files"
        echo "[WARN] DO NOT RE-RUN THIS SCRIPT"
        exit 1
      fi
      printf "\t[SUCCESS] Removing ${f}\n"
      rm ${f}
    done
  done
  cd -
  echo "Successful download & extraction"
  echo "Verifying DB"
  ${binary_dir}/blastdbcheck -dir ${out_path}

  echo "Export directory of db files as BLASTDB and run blast executables"
  echo "Blast executables: ${binary_dir}"
  echo ""

  echo "[Testing simple sequence]"
  test_file=".blast_test.fa"
  printf ">test\ntgcaccaaacatgtctaaagctggaaccaaaattactttctt\n" > ${test_file}
  test_out_file=".results.out"

  export_cmd="export BLASTDB=$(realpath ${out_path})"
  echo "${export_cmd}"
  eval ${export_cmd}

  blast_test_cmd="${binary_dir}/blastn -db ${db_name} -query test.fa -out ${test_out_file}"
  echo ${blast_test_cmd}
  eval ${blast_test_cmd}
  if [[ $? -eq 0 ]]; then
    echo "SUCCESS"
  else
    echo "FAIL"
  fi
  rm ${test_file} ${test_out_file}
  echo ""
else
  echo ""
  echo "For more info, see https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html"
  echo ""
  echo "Use this script to download databases"
  echo "perl ${download_script} <DB_NAME>"
  echo "  e.g."
  echo "    perl ${download_script} ref_euk_rep_genomes 	# Downloads all eukaryotic genomes"
fi
echo "Done."