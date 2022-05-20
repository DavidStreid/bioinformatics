#!/bin/bash
# Converts input BAM into paired-end FQ reads separated by chromosome

BAM=$1

base=$(basename ${BAM} | rev | cut -d'.' -f1- | rev)

mkdir -p ${base}
cd ${base}

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX; do
  base_x="${base}___${chr}"
  chrx_sam=${base_x}.bam
  chrx_sorted_sam=${base_x}.qsort.bam
  chrx_fq1=${base_x}.r1.fq
  chrx_fq2=${base_x}.r2.fq

  echo "Extracting ${chr} alignments into SAM=${chrx_sam}"
  samtools view -h ${BAM} "${chr}" > ${chrx_sam}
  wc -l ${chrx_sam}

  echo "Sorting ${chr} SAM=${chrx_sorted_sam}"
  samtools sort -n -o ${chrx_sorted_sam} ${chrx_sam}

  echo "Extracting ${chr} reads in FQ=${chrx_fq1} FQ2=${chrx_fq2}"
  bedtools bamtofastq -i ${chrx_sorted_sam} -fq ${chrx_fq1} -fq2 ${chrx_fq2}
done

cd -
