#!/bin/bash

RUN_FOLDER=$1 
SS=$2	

echo "SS=${SS} RUN_FOLDER=${RUN_FOLDER}"
if [[ ! -f ${SS} || ! -d ${RUN_FOLDER} ]]; then
  echo "[ERROR] Valid Run Folder (1) & Samplesheet (2) are needed" 
  exit 1
fi

RUN=$(basename ${RUN_FOLDER})
LOG="out_${RUN}.log"
ERR="out_${RUN}.err"
OUTPUT_DIR="$(pwd)/DEMUX_${RUN}"

echo "Running..."
nohup bcl-convert --bcl-input-directory ${RUN_FOLDER} \
  --output-directory ${OUTPUT_DIR} \
  --sample-sheet ${SS} > ${LOG} 2> ${ERR}
echo "Done."
