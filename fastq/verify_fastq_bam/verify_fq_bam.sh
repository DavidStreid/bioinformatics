#!/bin/bash

R1=$1
R2=$2
BAM=$3
OUT=$4

if [[ ! -f ${R1} || ! -f ${R2} || ! -f ${BAM} ]]; then
  echo "Invalid FQ/BAM"
  exit 1
fi

if [[ -z ${OUT} ]]; then
  OUT=stats.tsv
  printf "STATUS\tR1_READS\tR2_READS\tPAIRED_READS\tFQ_R1\tFQ_R2\tBAM\n" > ${OUT}
fi

for tool in samtools seqkit; do
  which ${tool}
  if [[ $? -ne 0 ]]; then
    echo "[ERROR] ${tool} isn't in path"
    exit 1
  fi
done

echo "#FASTQ"
echo "R1=${R1}"
echo "R2=${R2}"
r1_file=$(basename ${R1} | rev | cut -d'.' -f1- | rev)
r1_out="fq_stats___${r1_file}.tsv"
r2_file=$(basename ${R2} | rev | cut -d'.' -f1- | rev)
r2_out="fq_stats___${r2_file}.tsv"
seqkit stats ${R1} > ${r1_out}
seqkit stats ${R2} > ${r2_out}
r1_num_reads=$(grep "${r1_file}" "${r1_out}" | sed 's/   */\t/g' | cut -f4 | sed 's/,//g')
r2_num_reads=$(grep "${r2_file}" "${r2_out}" | sed 's/   */\t/g' | cut -f4 | sed 's/,//g')
total_fq_reads=$(( r1_num_reads + r2_num_reads ))

echo ""
b_file=$(basename ${BAM} | rev | cut -d'.' -f1- | rev)
b_stats="bam_stats___${b_file}.txt"
echo "#BAM"
echo "${BAM}"
samtools flagstat ${BAM} > ${b_stats}

num_paired_reads=$(grep "paired in sequencing" ${b_stats} | cut -d' ' -f1)

echo "STATS"
echo "R1_NUM_READS=${r1_num_reads}"
echo "R2_NUM_READS=${r2_num_reads}"
echo "BAM_PAIRED_READS=${num_paired_reads}"

status=""
if [[ ${num_paired_reads} -eq ${total_fq_reads} ]]; then
  status="LOSSLESS"
  echo "status=${status}"
  printf "\tFQ_READS=${total_fq_reads}\n"
  printf "\tBAM_PAIRED_READS=${num_paired_reads}\n"
else
  status="LOSSY"
  echo "[WARN] ${status} - Unpaired FASTQs"
  printf "\tFQ_READS=${total_fq_reads}\n"
  printf "\tBAM_PAIRED_READS=${num_paired_reads}\n"
fi

 echo "${status}\t${r1_num_reads}\t${r2_num_reads}\t${num_paired_reads}\t${R1}\t${R2}\t${BAM}\n" >> ${OUT}