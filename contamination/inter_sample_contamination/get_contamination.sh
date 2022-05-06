#!/bin/bash

BAM=$1 # Indexing required, 'samtools index ${BAM}'

if [[ ! -f ${BAM} ]]; then
  echo "Invalid BAM: '${BAM}'"
  exit 1
fi

sample="$(basename ${BAM} | cut -d'.' -f1)"

mkdir -p ${sample}
cd ${sample}
NEW_BAM=${sample}_temp.bam
if [[ ! -f ${NEW_BAM}.bai ]]; then
  # INDEX needs to be created 
  samtools sort ${BAM} > ${NEW_BAM}
  samtools index ${NEW_BAM}
fi

VCF="somatic-hg38_af-only-gnomad.hg38.vcf" 		# https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
REF="hg38.fa" 						                    # ftp://hgdownload.soe.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/bigZips/hg38.fa.gz
BED="hg38.bed"	                              # include all regions			# create


PILUEP_OUTPUT="${sample}_pileups___$(date +%s%N).table"

CMD="PATH/TO/gatk GetPileupSummaries "
CMD+="-I ${NEW_BAM} "
CMD+="-O ${PILUEP_OUTPUT} "
CMD+="-V ${VCF} "
CMD+="-L ${BED} "
CMD+="-R ${REF} "

P_LOG="log_${sample}_GetPileupSummaries___$(date +%s%N).out"
echo ${CMD}
echo ${CMD} > ${P_LOG}
eval ${CMD} >> ${P_LOG} 2>&1

if [[ -f ${NEW_BAM} ]]; then
  # Remove these temporary BAM files
  rm ${NEW_BAM}
fi

if [[ ! -f ${PILUEP_OUTPUT} ]]; then
  echo "[FAILED] pileups (${P_LOG})"
  echo ""
  exit 1
fi

CONTAMINATION_OUTPUT="${sample}_contamination___$(date +%s%N).table"
CMD="PATH/TO/gatk CalculateContamination "
CMD+="-I ${PILUEP_OUTPUT} "
CMD+="-O ${CONTAMINATION_OUTPUT}"

C_LOG="log_${sample}_CalculateContamination___$(date +%s%N).out"
echo ${CMD}
echo ${CMD} > ${C_LOG}
eval ${CMD} >> ${C_LOG} 2>&1

if [[ ! -f ${CONTAMINATION_OUTPUT} ]]; then
  echo "[FAILED] pileups (${C_LOG})"
  echo ""
  exit 1
fi

echo "[SUCCESS] (${CONTAMINATION_OUTPUT})"
echo ""

cd -
