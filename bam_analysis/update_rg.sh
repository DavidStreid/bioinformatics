#!/bin/bash

ID=$1
SAM=$2
SFX="bam"
if [[ -z ${ID} ]]; then
  echo "Need ID"
  exit 1
fi
if [[ ! -f ${SAM} ]]; then
  echo "Need SAM"
  exit 1
fi

FMT=$(echo ${SFX} | awk '{print toupper($0)}')
OUT="$(basename ${SAM} | rev | cut -f2- -d'.' | rev).${SFX}"
printf "FMT=${FMT}\tOUT=${OUT}\n"
samtools addreplacerg -w -r "ID:${ID}\tLB:Seq\tPL:Illumina\tSM:${ID}" ${SAM} -O ${FMT} -o ${OUT}
