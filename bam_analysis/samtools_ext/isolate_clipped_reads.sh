#!/bin/bash
# Isolates hard/soft clipped reads in SAM file

options="\tfilter - H=Hard\tS=Soft\tHS=Hard & Soft\n"
help_string="./isolate_clipped_reads.sh -b {BAM} [-f filter]\n${options}"

while getopts ":b:f:h" opt; do
    case $opt in
        b) INPUT_BAM=${OPTARG}
        ;;
        f) FILTER=${OPTARG}
        ;;
        h) printf "${help_string}"
        ;;
    esac 
done

if [[ ! -f ${INPUT_BAM} ]]; then
  echo "Please provide a valid BAM. INPUT_BAM=${INPUT_BAM}"
  printf "${help_string}"
  exit 1
fi

if [[ -z ${FILTER} ]]; then
  FILTER=HS # Default to both Hard & Soft-Clipped Reads
fi

REGEX="[0-9]+[${FILTER}].+[${FILTER}]$"

echo "BAM=${INPUT_BAM}"
echo "FILTER=${FILTER}"
printf "${options}"


OUTPUT_BAM="$(echo ${INPUT_BAM} | cut -d'.' -f1)_clipped.sam"

SC_RGID_FILE=soft_clipped_rgid.txt

# https://www.biostars.org/p/137461/#137631
echo "Extracting RGIDs of clipped BAM reads: ${SC_RGID_FILE}"
echo "REGEX=${REGEX}"
samtools view ${INPUT_BAM} | awk -v REGEX="${REGEX}" '$6~/'"$REGEX"'/{print $1}' | sort -k1,1 | uniq > ${SC_RGID_FILE}

echo "Found $(wc -l ${SC_RGID_FILE} | grep -oE "[0-9]+") clipped RGIDs"

echo "Writing ${OUTPUT_BAM}"
printf "\theaders...\n"
samtools view -H ${INPUT_BAM} > ${OUTPUT_BAM}

printf "\treads...\n"
samtools view ${INPUT_BAM} | fgrep -w -f ${SC_RGID_FILE} >> ${OUTPUT_BAM}

COLLATED_BAM="$(echo ${OUTPUT_BAM} | cut -d'.' -f1)___collated.bam"
echo "Writing Collated BAM=${COLLATED_BAM}"
samtools collate -o ${COLLATED_BAM} ${OUTPUT_BAM}

echo "Done."
