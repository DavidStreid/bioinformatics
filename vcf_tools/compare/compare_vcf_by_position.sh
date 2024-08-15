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

same_in_1=$(bedtools intersect -u -a  ${vcf1} -b ${vcf2} | wc -l)
same_in_2=$(bedtools intersect -u -a  ${vcf1} -b ${vcf2} | wc -l)
new_diffs=$(bedtools subtract -a ${vcf1} -b ${vcf2} | wc -l)
old_diffs=$(bedtools subtract -a ${vcf2} -b ${vcf1} | wc -l)
 

echo "vcf1_same=${same_in_1}"
echo "vcf2_same=${same_in_2}"
echo "vcf1_only=${new_diffs}"
echo "vcf2_only=${old_diffs}"
jaccard=$(bedtools jaccard -b $vcf1 -a $vcf2 | cut -f3)
echo "${jaccard}" | sed 's/ /=/g' # jaccard is [0,1] - 1 is exactly same positions, see https://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html
