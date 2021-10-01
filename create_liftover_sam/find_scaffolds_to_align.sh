#!/bin/bash

INPUT_REF=$1

LOCATION=$(dirname "$0")
SCAFFOLD_MAP_DIR="${LOCATION}/alternate_scaffold_mappings"

if [[ ! -f ${INPUT_REF} ]]; then
  echo "Please provide a valid reference file"
  exit 1
fi

num_contigs=$(cat ${INPUT_REF} | grep ">" | wc -l)
echo "Writing ${num_contigs} chr*.fa files for ${INPUT_REF}"
ref_to_scaffold_split_basename="chr"
faSplit sequence ${INPUT_REF} ${num_contigs} ${ref_to_scaffold_split_basename}
for f in ${ref_to_scaffold_split_basename}*.fa; do
  fname=$(head -1 ${f} | cut -d' ' -f1 | sed 's/>//')
  mv ${f} ${fname}.fa
done

echo "Grabbing all scaffolds from ${INPUT_REF}"
input_scaffolds=$(cat ${INPUT_REF} | grep ">" | cut -d' ' -f1 | sed 's/>//' | grep -E ".{3,}")
sMaps=$(find ${SCAFFOLD_MAP_DIR} -type f)
mapped_scaffolds_suffix="_map_scaffold.txt"

for map in ${sMaps}; do
  map_name=$(basename ${map} | cut -d'.' -f1)
  echo "Checking ${map_name}..."
  scaffolds_to_map_file="${map_name}${mapped_scaffolds_suffix}"
  touch ${scaffolds_to_map_file}
  for scaffold in ${input_scaffolds}; do
    grep "${scaffold}" ${map} >> ${scaffolds_to_map_file}
  done
  num_matches=$(wc -l ${scaffolds_to_map_file} | rev | cut -d' ' -f2 | rev)
  if [[ ${num_matches} -eq 0 ]]; then
    printf "\t... no matches\n"
    rm ${scaffolds_to_map_file}
  else
    printf "\t...${num_matches} matches\n"
    while IFS= read -r match; do
      qry=$(echo ${match} | cut -d' ' -f2)
      ref=$(echo ${match} | cut -d' ' -f6)
      echo "bwa mem -M -t 40 ${ref}.fa ${qry}.fa > q${qry}_r${ref}.sam"
    done < "${scaffolds_to_map_file}"
  fi
done
