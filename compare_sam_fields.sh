#!/bin/bash
# Compares the numerical fields of two SAM files created from the same set of FASTQ files
# Meant to compare aligners
# First intended use was to compare the MAPQ scores between a DRAGEN & BWA-MEM aligned BAM

B1=$1
B2=$2
IDX=$3 # MAPQ SCORE

if [[ -z ${B1} || -z ${B2} ]]; then
  echo "Please specify two .bam files and an optional index for a numerical column (e.g. 4: POSition, 5: MAPQ)"
  echo "e.g. compare_sam_fields.sh BAM1 BAM2 4"
  exit 1
fi

if [[ -z ${IDX} ]]; then
  IDX=5
  echo "Parsing out 5th column (MAPQ) of SAM files"
fi

B1_NAME="S1___$(echo $(basename ${B1}) | cut -d'.' -f1)"
B2_NAME="S2___$(echo $(basename ${B2}) | cut -d'.' -f1)"

B1_EXTRACT="EXTRACT___${B1_NAME}.csv"
B2_EXTRACT="EXTRACT___${B2_NAME}.csv"

echo "Extracting QNAME (1) & value at (${IDX}): ${B1_NAME} & ${B2_NAME}"
# Extract QNAME and additional SAM header field
samtools view ${B1} | cut -f1,2,${IDX} > ${B1_EXTRACT}
samtools view ${B2} | cut -f1,2,${IDX} > ${B2_EXTRACT}

echo "[CHECK 1] Same number of entries extracted"
B1_NUM=$(wc -l ${B1_EXTRACT} | cut -d' ' -f1)
B2_NUM=$(wc -l ${B2_EXTRACT} | cut -d' ' -f1)
if [[ ${B1_NUM} -ne ${B2_NUM} ]]; then
  echo "Different number of QNAME entires between input SAM files: ${B1_NUM} != ${B2_NUM}"
  exit 1
fi

echo "[CHECK 2] Same QNames in each SAM"
B1_QNAMES="QNAMES___${B1_NAME}.txt"
B2_QNAMES="QNAMES___${B2_NAME}.txt"
cat ${B1_EXTRACT} | cut -f1,2 | sort >> ${B1_QNAMES}
cat ${B2_EXTRACT} | cut -f1,2 | sort >> ${B2_QNAMES}
if [[ ! -z $(diff ${B1_QNAMES} ${B2_QNAMES}) ]]; then
  echo "Detected different QNAMES"
  diff ${B1_QNAMES} ${B2_QNAMES}
  exit 1
fi

RESULTS=RESULTS___${B1_NAME}_${B2_NAME}.csv
echo "Comparing selected values. Writing to ${RESULTS}"
while IFS= read -r b1_line
do
  read_id=$(echo ${b1_line} | cut -f1,2)
  b2_line=$(cat ${B2_EXTRACT} | grep "${read_id}")
  b1_val=$(echo ${b1_line} | cut -f3)
  b2_val=$(echo ${b2_line} | cut -f3)
  diff=$(expr ${b1_val} - ${b2_val})
  echo "${read_id},${b1_val},${b2_val},${diff}" >> ${RESULTS}
done < ${B1_EXTRACT}

echo "DONE"
