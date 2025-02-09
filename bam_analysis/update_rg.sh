#!/bin/bash
# Updates the readgroup of an input SAM to the input ID

ID=$1
SAM=$2
SFX=$3 # "bam" or "cram"
PL="Illumina"
if [[ -z ${ID} ]]; then
  echo "Need ID"
  exit 1
fi
if [[ ! -f ${SAM} ]]; then
  echo "Need SAM"
  exit 1
fi
if [[ -z ${SFX} ]]; then
  SFX=$(basename ${SAM} | rev | cut -f1 -d'.' | rev)
  printf "UPDATING SFX\t${SFX}\n"
fi

FMT=$(echo ${SFX} | awk '{print toupper($0)}')
OUT="$(basename ${SAM} | rev | cut -f2- -d'.' | rev).${SFX}"
printf "FMT=${FMT}\tOUT=${OUT}\tPL=${PL}\n"
samtools addreplacerg -@ 16 -w -r "ID:${ID}\tLB:Seq\tPL:${PL}\tSM:${ID}" ${SAM} -O ${FMT} -o ${OUT}
