#!/bin/bash

FROM_REF=$1
TO_REF=$2

if [[ ! -f ${FROM_REF} || ! -f ${TO_REF} ]]; then
  echo "Invalid Reference: ${FROM_REF}, ${TO_REF}"
  exit 1
fi

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
  echo ${INPUT_CMD}  >> ${LOG_FILE}
  eval ${INPUT_CMD} >> ${LOG_FILE} 2>&1
}

create_splits_dir() {
  ref_file=$1
  ref_dir=$(dirname ${ref_file})
  ref_base=$(basename ${ref_file} | cut -d'.' -f1)
  num_contigs=$(cat ${ref_file} | grep ">" | wc -l)
  log "Writing ${num_contigs} chr*.fa files for ${ref_file}"

  run_cmd "cd ${ref_dir}"
  ref_splits_dir="${WORK_DIR}/${ref_base}/splits"
  ref_size_dir="${WORK_DIR}/${ref_base}/size"
  ref_lft_dir="${WORK_DIR}/${ref_base}/lft"
  ref_2bit_dir="${WORK_DIR}/${ref_base}/2bit"
  scaffold_lft_dir="${ref_lft_dir}/scaffolds"

  mkdir -p ${ref_splits_dir} && \
    mkdir -p ${ref_lft_dir} && \
    mkdir -p ${ref_size_dir} && \
    mkdir -p ${ref_2bit_dir} && \
    mkdir -p ${scaffold_lft_dir}

  ref_lft_fname="lift_${ref_base}"
  run_cmd "faSplit sequence ${ref_file} ${num_contigs} -lift=${ref_lft_fname}"

  for chr in chr*.fa; do
    chr_base=$(basename ${chr} | cut -d'.' -f1)
    chr_2bit="${chr_base}.2bit"
    chr_size="${chr_base}.size"
    chr_lift="${chr_base}.lft"
    lift_prefix="lift_${chr_base}_"
    # We need to create the .lft files to change the coordinates after alignment
    log "Lift Creation: ${chr} -> ${chr_2bit} -> ${chr_size} -> ${chr_lift}"
    run_cmd "faToTwoBit ${chr} ${chr_2bit}"
    run_cmd "twoBitInfo ${chr_2bit} ${chr_size}"

    log "${chr_size}" 
    num_bases=$(cat ${chr_size} | cut -f2)
    log "${num_bases}" 
    run_cmd "faSplit size ${chr} ${num_bases} -lift=${chr_lift} ${lift_prefix}"
    num_files=$(ls -1 | grep "${lift_prefix}" | wc -l)
    if [[ ${num_files} -ne 1 ]]; then
      log "[ERROR] - Too many lift entries for scaffold: ${chr_base}"
      exit 1
    fi
    log "[SUCCESS] Created ${chr_lift}"
  done

  run_cmd "mv chr*.fa ${ref_splits_dir}"
  run_cmd "mv *.size ${ref_size_dir}"
  run_cmd "mv chr*.lft ${ref_lft_dir}"
  run_cmd "mv chr*.2bit ${ref_2bit_dir}"

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

blat_query() {
  to_splits_dir_param=$1
  from_2bit_param=$2
  log "to_splits_dir_param=${to_splits_dir_param} from_2bit_param=${from_2bit_param}"

  split_ref_files=$(find ${to_splits_dir_param} -type f -name "chr*.fa")
  for chr_fa in ${split_ref_files}; do
    log "Processing ${chr_fa}"
    chr_base=$(basename ${chr_fa} | cut -d'.' -f1)
    run_cmd "blat ${from_2bit_param} ${chr_fa} ${chr_base}.psl -tileSize=12 -minScore=100 -minIdentity=98" # -fastMap
  done
  ref_base_name=$(basename ${to_splits_dir_param})
  PSL_DIR=${WORK_DIR}/psl/${ref_base_name}
  mkdir -p ${PSL_DIR}
  log "Moving *.psl to ${PSL_DIR}"
   
  mv "chr*.psl" ${PSL_DIR}
  echo ${PSL_DIR}
}

echo "Splitting ${TO_REF}..."
to_splits_dir=$(create_splits_dir ${TO_REF})
printf "\t... wrote splits to ${to_splits_dir}\n"
echo "Creating 2bit of ${FROM_REF}..."
from_2bit=$(create_2bit_from_ref ${FROM_REF})
printf "\t... wrote 2bit to ${from_2bit}"

echo "BLAT query sequences of $(basename ${to_splits_dir}) against the 2bit..."
dir_blat_querys_against_FROM=$(blat_query ${to_splits_dir} ${from_2bit})
echo "... wrote .psl files to ${dir_blat_querys_against_FROM}"

# create_splits ${TO_REF}
# create_2bit_from_ref ${FROM_REF}
# twoBitInfo ${ref_2bit_name} chrom.sizes



# WORKS
# mkdir splits
# cd splits
# faSplit sequence ${REF} ${num_contigs} chr
# cd -

# Create 2bit genome sequence representation
# cd ${REF_DIR}
# REF_2BIT_NAME="${REF_BASE}.2bit"
# faToTwoBit ${REF} ${REF_2BIT_NAME}
# twoBitInfo ${REF_2BIT_NAME} chrom.sizes


