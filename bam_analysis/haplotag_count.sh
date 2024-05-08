#!/bin/bash
# Quick check to see how many reads contain haplotype (HP) identifiers

SAM=$1

echo "SAM=${SAM}"
if [[ ! -f ${SAM} ]]; then
  echo "Invalid File"
  exit 1
fi

log=$(basename ${SAM}).out
touch ${log}
samtools view ${SAM} | grep -oE "HP:i:[12]" | wc -l >> ${log}
echo "  Finished HP count"
samtools view ${SAM} | wc -l >> ${log}
echo "  Finished total count"
cat ${log}
