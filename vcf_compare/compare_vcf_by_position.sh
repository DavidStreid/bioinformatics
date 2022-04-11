#!/bin/bash
# Compares two VCF files, treating their position as IDs

vcf1=$1
vcf2=$2

echo "vcf1=${vcf1}"
echo "vcf2=${vcf2}"

if [[ ! -f ${vcf1} || ! -f ${vcf2} ]]; then
  echo "Two valid files are required. Exiting."
  exit 1
fi

same=$(bedtools intersect -u -a  ${vcf1} -b ${vcf2} | wc -l)
new_diffs=$(bedtools subtract -a ${vcf1} -b ${vcf2} | wc -l)
old_diffs=$(bedtools subtract -a ${vcf2} -b ${vcf1} | wc -l)

echo "same=${same}"
echo "new=${new_diffs}"
echo "old=${old_diffs}"
