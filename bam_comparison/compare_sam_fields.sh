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
samtools view ${B1} | cut -f1,2,${IDX} | sort > ${B1_EXTRACT}
samtools view ${B2} | cut -f1,2,${IDX} | sort > ${B2_EXTRACT}

echo "[CHECK 1] Same number of entries extracted"
B1_NUM=$(wc -l ${B1_EXTRACT} | cut -d' ' -f1)
B2_NUM=$(wc -l ${B2_EXTRACT} | cut -d' ' -f1)
if [[ ${B1_NUM} -ne ${B2_NUM} ]]; then
  printf "\t[WARNING]Different number of QNAME entires between input SAM files: ${B1_NUM} != ${B2_NUM}\n"
fi

echo "[CHECK 2] Same QNames in each SAM"
B1_QNAMES="QNAMES___${B1_NAME}.txt"
B2_QNAMES="QNAMES___${B2_NAME}.txt"
cat ${B1_EXTRACT} | cut -f1,2 | sort >> ${B1_QNAMES}
cat ${B2_EXTRACT} | cut -f1,2 | sort >> ${B2_QNAMES}
if [[ ! -z $(diff ${B1_QNAMES} ${B2_QNAMES}) ]]; then
  printf "\t[WARNING] Detected different QNAMES\n"
fi

RESULTS=RESULTS___${B1_NAME}_${B2_NAME}.csv
echo "QNAME,FLAG_1,FLAG_2,V1,V2,DIFF" >> ${RESULTS}

echo "Comparing selected values. Writing to ${RESULTS}"
while IFS= read -r b1_line
do
  echo "${b1_line}"
  qname_1=$(echo ${b1_line} | cut -d' ' -f1)
  flag_1=$(echo ${b1_line} | cut -d' ' -f2)

  r1=$(cat ${B2_EXTRACT} | grep "${qname_1}" | head -1)
  r2=$(cat ${B2_EXTRACT} | grep "${qname_1}" | tail -1)

  r1_flag=$(echo $r1 | cut -d' ' -f2)
  r2_flag=$(echo $r2 | cut -d' ' -f2)

  printf "\tChoosing read for flag_1: ${flag_1}\n"
  printf "\tr1_flag: ${r1_flag}\n"
  printf "\tr2_flag: ${r2_flag}\n"

  if [[ "${flag_1}" -gt "127" ]]; then
    # input flag indicates input read is R1 
    if [[ "${r1_flag}" -gt "127" ]]; then
      printf "\tChoosing r1 ($r1_flag): $r1\n"
      b2_line=$r1
      flag_2=$r1_flag
    else
      printf "\tChoosing r2 ($r2_flag): $r2\n"
      b2_line=$r2
      flag_2=$r2_flag
    fi
  else
    # input flag indicates input read is R2
    if [[ "${r1_flag}" -lt "128" ]]; then
      printf "\tChoosing r1 ($r1_flag): $r1\n"
      b2_line=$r1
      flag_2=$r1_flag
    else
      printf "\tChoosing r2 ($r2_flag): $r2\n"
      b2_line=$r2
      flag_2=$r2_flag
    fi    
  fi

  if [[ -z ${b2_line} ]]; then
    printf "\t[SKIPPING] Could not find ${read_id} in ${B2_EXTRACT}"
    continue
  else
    printf "\tComparing the following two lines\n"
    printf "\t${B1_EXTRACT}\n"
    printf "\t${b1_line}\n"
    printf "\t${B2_EXTRACT}\n"
    printf "\t${b2_line}\n"
  fi
  b1_val=$(echo ${b1_line} | cut -d' ' -f3)
  b2_val=$(echo ${b2_line} | cut -d' ' -f3)
  diff=$(expr ${b1_val} - ${b2_val})
  results_line="${qname_1},${flag_1},${flag_2},${b1_val},${b2_val},${diff}"
  printf "\t${results_line}\n"
  echo ""
  echo "${results_line}" >> ${RESULTS}
done < ${B1_EXTRACT}

echo "DONE"
