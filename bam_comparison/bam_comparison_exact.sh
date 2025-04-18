#!/bin/bash

B1=$1
B2=$2

b1_base=$(basename ${B1} | sed 's/.bam$//g')
b2_base=$(basename ${B2} | sed 's/.bam$//g')
for chrom in $(samtools idxstats ${B1} | grep -Ev "random|chrUn|\*" | cut -f1); do
  echo ${chrom}
  c1=${chrom}.${b1_base}.sam
  c2=${chrom}.${b2_base}.sam
  samtools view -h ${B1} "${chrom}" | sort > ${chrom}.${b1_base}.sam
  samtools view -h ${B2} "${chrom}" | sort > ${chrom}.${b2_base}.sam
  diff ${c1} ${c2} > diff_${chromt}.txt
done
