#!/bin/bash

f1=s1.vcf
f2=s2.vcf

if [[ ! -f ${f1} || ! -f ${f2} ]]; then
  echo "Need f1 & f2"
  exit 1
fi
info=$(grep -A 99999999999 "CHROM" ${f1} ${f2}  | cut -f8 | sed 's/;/\n/g' | cut -d'=' -f1 | sort | uniq)
format=$(grep -A 99999999999 "CHROM" ${f1} ${f2}  | cut -f9 | sed 's/:/\n/g' | sort | uniq | grep -Ev FORMAT | grep -Ev "^$")

echo "##fileformat=VCFv4.1"
echo "##source=VCF_MERGE"
for x in ${info}; do
  grep "<ID=${x}," ${f1} ${f2} | cut -d':' -f2- | sort | uniq  |  grep INFO
done
for x in ${format}; do
  grep  "<ID=${x}," ${f1} ${f2} | cut -d':' -f2- | sort | uniq  | grep FORMAT
done
