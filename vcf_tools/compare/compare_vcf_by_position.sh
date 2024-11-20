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

shared_file="shared.tsv"
vcf1_only="vcf1_only.tsv"
vcf2_only="vcf2_only.tsv"
cmds=("bedtools intersect -f 1 -r -a ${vcf1} -b ${vcf2} > ${shared_file}"   # same exact position
"bedtools subtract -f 1 -a ${vcf1} -b ${vcf2} > ${vcf1_only}"               # VCF1 variants with NO overlap
"bedtools subtract -f 1 -a ${vcf2} -b ${vcf1} > ${vcf2_only}")              # VCF2 variants with NO overlap
for cmd in "${cmds[@]}"; do
  echo "${cmd}"
  eval "${cmd}"
done


shared_ct=$(cat ${shared_file} | wc -l)
old_diffs=$(cat ${vcf1_only} | wc -l)
new_diffs=$(cat ${vcf2_only} | wc -l)
 

echo "shared=${shared_ct}"
echo "vcf1_only=${new_diffs}"
echo "vcf2_only=${old_diffs}"
jaccard=$(bedtools jaccard -b $vcf1 -a $vcf2 | cut -f3)
echo "${jaccard}" | sed 's/ /=/g' # jaccard is [0,1] - 1 is exactly same positions, see https://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html
