#!/bin/bash

print_usage() {
  echo "./downsample.sh -f [input_bam] -r [num_reads]"
  printf "\t./downsample.sh -f sample.fastq -r 1000\n"
}

while getopts 'f:r:' flag; do
  case "${flag}" in
    f) INPUT_FASTQ="${OPTARG}" ;;     # Reference to create liftover file for, e.g. GRCh37
    r) NUM_READS="${OPTARG}" ;;          # Option to submit the blat job to lsf or do it locally
    *) print_usage
       exit 1 ;;
  esac
done

if [[ -z ${INPUT_FASTQ} || ! -f ${INPUT_FASTQ} ]]; then
  echo "[ERROR] Invalid FASTQ: ${INPUT_FASTQ}"
  print_usage
  exit 1
fi

if [[ -z ${NUM_READS} ]]; then
  echo "[ERROR] Please specify a number of reads"
  print_usage
  exit 1
fi

FASTQ_FILE=$(basename ${INPUT_FASTQ})
FASTQ_BASE=$(echo ${FASTQ_FILE} | cut -d'.' -f1)
DOWNSAMPLED_BAM="DOWNSAMPLED_${FASTQ_BASE}.fastq"

# -2 uses "two-pass" mode, which consumes much less memory
#   Ref: https://github.com/lh3/seqtk/issues/140
seqtk sample -2 -s100 \
  ${INPUT_FASTQ} \
  ${NUM_READS} >> ${DOWNSAMPLED_BAM}


