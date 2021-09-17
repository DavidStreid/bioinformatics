#!/bin/bash

while getopts 's:t:l' flag; do
  case "${flag}" in
    s) FROM_REF="${OPTARG}" ;;
    t) TO_REF="${OPTARG}" ;;
    l) USE_LOCAL='true' ;;
    *) print_usage
       exit 1 ;;
  esac
done

if [[ ! -f ${FROM_REF} || ! -f ${TO_REF} ]]; then
  echo "Invalid Reference: ${FROM_REF}, ${TO_REF}"
  exit 1
fi

echo "Creating Liftover: ${FROM_REF} -> ${TO_REF} (USE_LOCAL=${USE_LOCAL})"
FROM_REF=$(realpath ${FROM_REF})
TO_REF=$(realpath ${TO_REF})

WORK_DIR=$(pwd)/temp
LOG_FILE=${WORK_DIR}/liftover.out
mkdir -p ${WORK_DIR}

log() {
  INPUT_CMD=$@
  echo ${INPUT_CMD}  >> ${LOG_FILE}
}

#########################################
# Executes and logs command
# Arguments:
#   INPUT_CMD - string of command to run, e.g. "picard CollectAlignmentSummaryMetrics ..."
#########################################
run_cmd () {
  INPUT_CMD=$@
  echo ${INPUT_CMD} >> ${LOG_FILE}
  eval "${INPUT_CMD}" >> ${LOG_FILE} 2>&1
}



create_setup_files() {
  ref_file=$1
  ref_2bit=$2
  setup_info_file=$3 # "setup_info.txt"
  ref_dir=$(dirname ${ref_file})
  ref_base=$(basename ${ref_file} | cut -d'.' -f1)
  num_contigs=$(cat ${ref_file} | grep ">" | wc -l)
  log "Writing ${num_contigs} chr*.fa files for ${ref_file}"

  run_cmd "cd ${ref_dir}"
  ref_splits_dir="${WORK_DIR}/${ref_base}/splits"
  ref_size_dir="${WORK_DIR}/${ref_base}/size"
  ref_lft_dir="${WORK_DIR}/${ref_base}/lft"
  ref_2bit_dir="${WORK_DIR}/${ref_base}/2bit"
  ref_psl_dir="${WORK_DIR}/${ref_base}/psl"

  mkdir -p ${ref_splits_dir} && \
    mkdir -p ${ref_lft_dir} && \
    mkdir -p ${ref_size_dir} && \
    mkdir -p ${ref_2bit_dir} && \
    mkdir -p ${ref_psl_dir}

  ref_to_scaffold_split_basename="chr"
  run_cmd "faSplit sequence ${ref_file} ${num_contigs} -lift=${ref_to_scaffold_split_basename}.lft ${ref_to_scaffold_split_basename}"

  JOB_ID_LIST_FILE="blat_lsf_jobs.txt"
  for chr in ${ref_to_scaffold_split_basename}*; do
    chr_base=$(basename ${chr} | cut -d'.' -f1)
    chr_2bit="${chr_base}.2bit"
    chr_size="${chr_base}.size"
    chr_lift="${chr_base}.lft"
    chr_blat="${chr_base}.psl"

    lift_prefix="lift_${chr_base}_"

    # We need to create the .lft files to change the coordinates after alignment
    log "Setup Files: ${chr} -> ${chr_2bit} -> ${chr_size} -> ${chr_lift} -> ${chr_blat}"
    run_cmd "faToTwoBit ${chr} ${chr_2bit}"
    run_cmd "twoBitInfo ${chr_2bit} ${chr_size}"

    num_bases=$(cat ${chr_size} | cut -f2)
    run_cmd "faSplit size ${chr} ${num_bases} -lift=${chr_lift} ${lift_prefix}"
    num_files=$(ls -1 | grep "${lift_prefix}" | wc -l)
    if [[ ${num_files} -ne 1 ]]; then
      log "[ERROR] - Too many lift entries for scaffold: ${chr_base}"
      exit 1
    fi

    chr_base=$(basename ${chr} | cut -d'.' -f1)
    BLAT_CMD="blat ${ref_2bit} ${chr} ${chr_blat} -tileSize=12 -minScore=100 -minIdentity=98"
    if [[ ${USE_LOCAL} == 'true' ]]; then
      run_cmd ${BLAT_CMD}
    else
      JOB_NAME="BLAT_${chr_base}"
      SUBMIT=$(bsub -J ${JOB_NAME} -o ${JOB_NAME} -n 10 -M 12 "${BLAT_CMD}")
      log ${SUBMIT}
      JOB_ID=$(echo $SUBMIT | egrep -o '[0-9]{5,}') # Parses out job id from output
      echo ${JOB_ID} >> ${JOB_ID_LIST_FILE}         # Save job id to wait on later
    fi

    log "[SUCCESS] Created ${chr_lift}"
  done

  if [[ -f ${JOB_ID_LIST_FILE} && USE_LOCAL == "false" ]]; then
    # LSF was chosen
    ALL_JOBS=$(cat ${JOB_ID_LIST_FILE} | tr '\n' ' ')
    log "Submitted Jobs (${EXECUTOR}): %s"
    log "${ALL_JOBS}"
    for job_id in ${ALL_JOBS}; do
      log "Waiting for ${job_id} to finish"
      bwait -w "ended(${job_id})" &
    done
  fi

  # Organize all outputs
  run_cmd "mv 'chr*.psl' ${ref_psl_dir}"
  run_cmd "mv chr*.fa ${ref_splits_dir}"
  run_cmd "mv *.size ${ref_size_dir}"
  run_cmd "mv chr*.lft ${ref_lft_dir}"
  run_cmd "mv chr*.2bit ${ref_2bit_dir}"

  # Write all files for each scaffold to a line in a file
  for split_file in ${ref_splits_dir}/${ref_to_scaffold_split_basename}*.fa; do
    bname=$(basename ${chr} | cut -d'.' -f1)
    bname_size=${ref_size_dir}/${bname}.size
    bname_lft=${ref_lft_dir}/${bname}.lft
    bname_2bit=${ref_2bit_dir}/${bname}.2bit
    bname_psl=${ref_psl_dir}/${bname}.psl

    log "Searching for ${split_file} ${bname_size} ${bname_lft} ${bname_2bit} ${bname_psl}"
    if [[ -f ${split_file} && -f ${bname_size} && -f ${bname_lft} && -f ${bname_2bit} && -f ${bname_psl} ]]; then
      echo "${bname} ${split_file} ${bname_size} ${bname_lft} ${bname_2bit} ${bname_psl}" >> ${setup_info_file}
    else
      echo "${bname} ERROR" >> ${setup_info_file}
      log "[ERROR] Setup of ${bname} "
    fi
  done

  cd - > /dev/null 	# Don't want to echo this out
  echo $(realpath ${ref_splits_dir})
}

create_2bit_from_ref() {
  ref_file=$1
  ref_dir=$(dirname ${ref_file})
  ref_base=$(basename ${ref_file} | cut -d'.' -f1)
  ref_2bit_name="${ref_base}.2bit"
  dir_for_2bit="${WORK_DIR}/2bit"
  mkdir -p ${dir_for_2bit}
  cd ${dir_for_2bit}
  run_cmd "faToTwoBit ${ref_file} ${ref_2bit_name}"
  cd - > /dev/null      # Don't want to echo this out
  echo "${dir_for_2bit}/${ref_2bit_name}"
}

echo "Creating 2bit of ${FROM_REF}..."
from_2bit=$(create_2bit_from_ref ${FROM_REF})
printf "\t... wrote 2bit to ${from_2bit}\n"

TO_SETUP_FILE="TO_SETUP_FILE.txt"
echo "Setup for ${TO_REF}. Writing to ${TO_SETUP_FILE}..."
create_setup_files ${TO_REF} ${from_2bit} ${TO_SETUP_FILE}
cat ${TO_SETUP_FILE} | grep ERROR
