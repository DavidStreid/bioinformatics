#!/bin/bash
# Wrapper to create the blastdb files from a list of taxonomic IDs IF NCBI has a mapping for that taxonomic ID

taxid_list=$@
sources_file="blastdb_sources.tsv"
log_file="log_blastdb_sources.out"

if [[ -z ${taxid_list} ]]; then
  echo "At least one TAXONOMIC ID is required. Exiting"
  exit 1
fi

echo "Locating sources for list of taxonomic IDs"
for taxid in ${taxid_list}; do
  echo "Finding source for ${taxid}"
  ./staxis_to_ncbi_assemlby.sh ${taxid} > ${sources_file}
  tail -1 ${sources_file}
done

echo "Creating BLASTDB files"
./create_blastdb_files.sh ${sources_file}