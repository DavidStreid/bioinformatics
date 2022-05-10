#!/bin/bash
# Script to find the best NCBI genome link for the input taxonomic ID
# Note: Info on columns of @TAXID_MAP_FILE at ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt

# Allowed DB files are viruses and human genome
DB_FILES=""
TAXID_MAP_FILE=assembly_summary_genbank.txt
TAXID="$1"

get_assembly_level() {
  # Returns the best ftp link for the input @type & @taxid. Selection Criteria below,
  #   - Latest version (11, "version_status")
  #   - assembly_level (12, highest level of assembly) - priority (in-order) for "Complete genome", "Chromosome", "Scaffold"
  #   - No exclusion reason (21, excluded_from_refseq) - this should be empty
  #   - Latest valid entry (get last)
  # @param type, refseq_category: e.g. "reference genome", "representative genome", "na"
  # @global TAXID, NCBI taxonomy identifier
  type=$@
  assembly="Complete genome"
  complete_available=$(awk -v TYPE="${type}" -v TAXID="${TAXID}" -v ASSEMBLY="${assembly}" \
    'BEGIN{FS="\t"} { \
      if($5==TYPE && $6==TAXID && $11=="latest" && $12==ASSEMBLY && $21=="") { print $20 } }' ${TAXID_MAP_FILE} \
        | tail -1 \
        | awk 'BEGIN{OFS=FS="/"}{print $0,$NF"_genomic.fna.gz"}')
  if [[ ! -z ${complete_available} ]]; then
    printf "${assembly}\t${complete_available}\n"
  else
    assembly="Chromosome"
    chr_available=$(awk -v TYPE="${type}" -v TAXID="${TAXID}"  -v ASSEMBLY="${assembly}" \
      'BEGIN{FS="\t"} { \
        if($5==TYPE && $6==TAXID && $11=="latest" && $12==ASSEMBLY && $21=="") { print $20 } }' ${TAXID_MAP_FILE} \
        | tail -1 \
        | awk 'BEGIN{OFS=FS="/"}{print $0,$NF"_genomic.fna.gz"}')
    if [[ ! -z ${chr_available} ]]; then
      printf "${assembly}\t${chr_available}\n"
    else
      assembly="Scaffold"
      scfld_available=$(awk -v TYPE="${type}" -v TAXID="${TAXID}" -v ASSEMBLY="${assembly}" \
        'BEGIN{FS="\t"} { \
          if($5==TYPE && $6==TAXID && $11=="latest" && $12==ASSEMBLY && $21=="") { print $20 } }' ${TAXID_MAP_FILE} \
          | tail -1 \
          | awk 'BEGIN{OFS=FS="/"}{print $0,$NF"_genomic.fna.gz"}')
      if [[ ! -z ${scfld_available} ]]; then
        printf "${assembly}\t${scfld_available}\n"
      fi
    fi
  fi
}

# Outline borrowed from this biostars post - https://www.biostars.org/p/306380/#306654
if [[ -z "$TAXID" ]]; then
  echo "A taxonomic ID is required to run. Usage is below -"
  printf "\t./staxis_to_ncbi_assemlby.sh \${TAX_ID}\n"
  exit 1
else
  # Check if tax_id is present in taxonomic databases
  if [[ ! -z ${DB_FILES} ]]; then
    db_file=$(grep -lE "${TAXID}\t" ${DB_FILES})
  fi
  if [[ ! -z ${db_file} ]]; then
    src=$(basename ${db_file})
    ftp="BLAST_DB\tn/a"
  else
    if [[ ! -e assembly_summary_genbank.txt ]]; then
      >&2 echo "Downloading ${TAXID_MAP_FILE}"
      wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/${TAXID_MAP_FILE}
    fi
    >&2 echo "identifying organism for TAXID=${TAXID}"
    ORGANISM=$(cut -f6,8 ${TAXID_MAP_FILE} | grep -E "^${TAXID}\t" | cut -f2 | head -1)
    >&2 echo "ORGANISM=${ORGANISM}"
    src="reference genome"
    >&2 echo "Checking for valid downloads where refseq_category=${src}"
    ftp=$(get_assembly_level ${src})
    if [[ -z ${ftp} ]]; then
      src="representative genome"
      >&2 echo "Checking for valid downloads where refseq_category=${src}"
      ftp=$(get_assembly_level ${src})
      if [[ -z ${ftp} ]]; then
        src="na"
        >&2 echo "Checking for valid downloads where refseq_category=${src}"
        ftp=$(get_assembly_level ${src})
        if [[ -z ${ftp} ]]; then
          src="NONE"
          ftp="NONE\tNONE"
        fi
      fi
    fi
  fi
  printf "${TAXID}\t${ORGANISM}\t${src}\t${ftp}\n"
fi
