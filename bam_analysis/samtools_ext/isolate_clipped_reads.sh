#!/bin/bash
# Isolates hard/soft clipped reads in SAM file

options="\t-f, file: SAM/BAM file\n\t-c, clipping: H=Hard\tS=Soft\tHS=Hard & Soft\n\t-p, paried-clipping: Flag to only include reads where both reads have the \"-c\" clipping"
help_string="./isolate_clipped_reads.sh -f file [-c clipping] [-p paired-clipping]\n${options}"

while getopts ":f:c:ph" opt; do
    case $opt in
        f) INPUT_BAM=${OPTARG}
        ;;
        c) CLIPPING_FILTER=${OPTARG}
        ;;
        p) PAIRED_CLIPPING_FILTER="TRUE"
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

# If filtering on both paired reads, check each RGID for two entries in our file of RGIDs w/ clipped-reads
tmp_file=".${CLIPPED_RGID_FILE}"
touch ${tmp_file}
if [[ ! -z ${PAIRED_CLIPPING_FILTER} ]]; then
  echo "Including only RGIDs w/ clipping of both paired-reads"
  for rgid in $(cat ${CLIPPED_RGID_FILE} | uniq); do
    if [[ 1 -ne $(grep "${rgid}" ${CLIPPED_RGID_FILE} | wc -l) ]]; then
      echo ${rgid} >> ${tmp_file}
    fi
  done
  num_filtered=$(wc -l ${tmp_file} | grep -oE "[0-9]+")
  echo "Found ${num_filtered} RGIDs with clipping of both pairs"
  if [[ 0 -eq ${num_filtered} ]]; then
    echo "No reads match the provided filters. Not writing BAM & exiting."
    exit 0
  fi 
else
  echo "Including all RGIDs w/ clipping of at least one paired-read"
  sort ${CLIPPED_RGID_FILE} | uniq > ${tmp_file}
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

echo "Done."

