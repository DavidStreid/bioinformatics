#!/bin/bash
# Compares two VCF files, treating their position as IDs
# NOTE - writes the output to a file

# TODO
# PATH_TO_BEDTOOLS=""
# PATH="${PATH}:${PATH_TO_BEDTOOLS}"

vcf1=$1
vcf2=$2

echo "vcf1=${vcf1}"
echo "vcf2=${vcf2}"

if [[ ! -f ${vcf1} || ! -f ${vcf2} ]]; then
  echo "Two valid files are required. Exiting."
  exit 1
fi


f_same_1="shared_vcf1.tsv"
f_same_2="shared_vcf2.tsv"
f_diff_1="vcf1_only.tsv"
f_diff_2="vcf2_only.tsv"
cmds=("bedtools intersect -u -a ${vcf1} -b ${vcf2} > ${f_same_1}"
"bedtools intersect -u -a ${vcf1} -b ${vcf2} > ${f_same_2}"
"bedtools subtract -a ${vcf1} -b ${vcf2} > ${f_diff_1}"
"bedtools subtract -a ${vcf2} -b ${vcf1} > ${f_diff_2}")
for cmd in "${cmds[@]}"; do
  echo "${cmd}"
  eval "${cmd}"
done


same_in_1=$(cat ${f_same_1} | wc -l)
same_in_2=$(cat ${f_same_2} | wc -l)
new_diffs=$(cat ${f_diff_1} | wc -l)
old_diffs=$(cat ${f_diff_2} | wc -l)
 

echo "vcf1_same=${same_in_1}"
echo "vcf2_same=${same_in_2}"
echo "vcf1_only=${new_diffs}"
echo "vcf2_only=${old_diffs}"
jaccard=$(bedtools jaccard -b $vcf1 -a $vcf2 | cut -f3)
echo "${jaccard}" | sed 's/ /=/g' # jaccard is [0,1] - 1 is exactly same positions, see https://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html
