#!/bin/bash
# Isolates hard/soft clipped reads in SAM file
# TODO - Add single-clipping option

options="\t-f, file: SAM/BAM file\n\t-c, clipping: H=Hard\tS=Soft\tHS=Hard & Soft\n\t-n, number of rgids: Max number of RGIDs to put in filtered BAM\n\t-p, paried-clipping: Flag to only include reads where both reads have the\n"
help_string="./isolate_clipped_reads.sh -f file [-c clipping] [-n num_reads] [-p] [-q]\n${options}"

while getopts ":f:c:n:pqh" opt; do
    case $opt in
        f) INPUT_BAM=${OPTARG}
        ;;
        c) CLIPPING_FILTER=${OPTARG}
        ;;
        n) NUM_READS=${OPTARG}
        ;;
        p) PAIRED_CLIPPING_FILTER="TRUE"
        ;;
        q) OUTPUT_FQ="TRUE"
        ;;
        h) printf "${help_string}" && exit 0
        ;;
    esac 
done

if [[ ! -f ${INPUT_BAM} ]]; then
  echo "Please provide a valid BAM. INPUT_BAM=${INPUT_BAM}"
  printf "${help_string}"
  exit 1
fi
if [[ -z ${CLIPPING_FILTER} ]]; then
  CLIPPING_FILTER=HS # Default to both Hard & Soft-Clipped Reads
fi
if [[ ! -z ${NUM_READS} ]]; then
  re='^[0-9]+$'
  if ! [[ ${NUM_READS} =~ $re ]] ; then
     echo "NUM_READS=${NUM_READS} is not a number. Exiting" >&2; exit 1
  fi
fi

echo "[INPUTS]"
echo "BAM=${INPUT_BAM}"
echo "CLIPPING_FILTER=${CLIPPING_FILTER}"
echo "PAIRED_CLIPPING_FILTER=${PAIRED_CLIPPING_FILTER}"
echo ""

# REF: https://www.biostars.org/p/137461/#137631
CLIPPED_RGID_FILE="${CLIPPING_FILTER}_clipped_rgid.txt"
echo "Extracting RGIDs of clipped BAM reads: ${CLIPPED_RGID_FILE}"
REGEX="[0-9]+[${CLIPPING_FILTER}].+[${CLIPPING_FILTER}]$"
echo "REGEX=\"${REGEX}\""
samtools view ${INPUT_BAM} | awk -v REGEX="${REGEX}" '$6~/'"$REGEX"'/{print $1}' | sort -k1,1 > ${CLIPPED_RGID_FILE}
echo "Found $(wc -l ${CLIPPED_RGID_FILE} | grep -oE "[0-9]+") clipped RGIDs"
echo ""

# If no limit was specified, take all reads
if [[ -z ${NUM_READS} ]]; then
  NUM_READS=$(wc -l ${CLIPPED_RGID_FILE})
fi

# If filtering on both paired reads, check each RGID for two entries in our file of RGIDs w/ clipped-reads
tmp_file=".${CLIPPED_RGID_FILE}"
touch ${tmp_file}
if [[ ! -z ${PAIRED_CLIPPING_FILTER} ]]; then
  added_reads=0
  echo "Including only RGIDs w/ clipping of both paired-reads"
  for rgid in $(cat ${CLIPPED_RGID_FILE} | uniq); do
    if [[ 1 -ne $(grep "${rgid}" ${CLIPPED_RGID_FILE} | wc -l) ]]; then
      echo ${rgid} >> ${tmp_file}
      added_reads=$((added_reads + 1))
      if [[ ${added_reads} -eq ${NUM_READS} ]]; then
        break
      fi
    fi
  done
  num_filtered=$(wc -l ${tmp_file} | grep -oE "[0-9]+")
  printf "\tAdded ${num_filtered} RGIDs with clipping of both pairs\n"
  if [[ 0 -eq ${num_filtered} ]]; then
    echo "\tNo reads match the provided filters. Not writing BAM & exiting.\n"
    exit 0
  fi 
else
  echo "Including all RGIDs w/ clipping of at least one paired-read"
  head -${NUM_READS} ${CLIPPED_RGID_FILE} | sort | uniq > ${tmp_file}
fi
mv ${tmp_file} ${CLIPPED_RGID_FILE}
echo ""

CLIPPED_READ_BAM="$(basename ${INPUT_BAM} | cut -d'.' -f1)_clipped.sam"
echo "Writing ${CLIPPED_READ_BAM}"
printf "\theaders...\n"
samtools view -H ${INPUT_BAM} > ${CLIPPED_READ_BAM}
printf "\treads...\n"
samtools view ${INPUT_BAM} | fgrep -w -f ${CLIPPED_RGID_FILE} >> ${CLIPPED_READ_BAM}
echo ""

COLLATED_BAM="$(echo ${CLIPPED_READ_BAM} | cut -d'.' -f1)_collated.bam"
echo "Shuffling/Grouping reads into collated BAM=${COLLATED_BAM}"
samtools collate -o ${COLLATED_BAM} ${CLIPPED_READ_BAM}
echo ""

if [[ ! -z ${OUTPUT_FQ} ]]; then
  FQ_FILE="$(echo ${COLLATED_BAM} | cut -d'.' -f1).fq"
  rm -f ${FQ_FILE}
  echo "Outputing FASTQ file: ${FQ_FILE}"
  samtools view ${COLLATED_BAM} | \
    while IFS= read line; do
      CIGAR=$(echo "${line}" | cut -f6)
      QNAME=$(echo "${line}" | cut -f1)
      ID="@${QNAME}___${CIGAR}"
      echo ${ID} >> ${FQ_FILE}
      echo "${line}" | cut -f10 >> ${FQ_FILE}
      echo "+" >> ${FQ_FILE}
      echo "${line}" | cut -f11 >> ${FQ_FILE}
    done
fi

echo "Done."

