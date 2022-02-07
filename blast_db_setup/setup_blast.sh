#!/bin/bash

version=$1
os_input=$2
db_name=$3

if [[ -z ${version} || ${version} == "l" ]]; then
  download_version=$(curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/ 2> /dev/null | grep LATEST | sed 's/.*-> //g')
else
  download_version=${version}
fi

if [[ -z ${os_input} ]]; then
  os="x64-linux"
else
  os=${os_input}
fi

untarred_file="ncbi-blast-${download_version}+"
tar_file="${untarred_file}-${os}.tar.gz"
echo "blast+ Version: ${download_version}"
echo "TAR file: ${tar_file}"
echo "os=${os_input}"
echo ""

if [[ -f ${tar_file} ]]; then
  echo "Detected downloaded file... skipping download"
elif [[ -f ${untarred_file} ]]; then
  echo "Detected extracted file... skipping download"
else
  CMD="curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${download_version}/${tar_file} -o ${tar_file}"
  echo ${CMD}
  eval ${CMD}
  if [[ $? -ne 0 ]]; then
    printf "\nDownload failed. Exiting...\n"
    exit 1
  else
    echo "\nDownload complete.\n"
  fi
fi

 
if [[ -d ${untarred_file} ]]; then
  echo "Nothing to untar - Using ${untarred_file}"
elif [[ -f ${tar_file} ]]; then
  echo "Extracting contents..."
  tar -zxvf ${tar_file}
  echo "Extraction done."
else
  echo "No tarred/untarred file. Exiting..."
  exit 1
fi

download_script=$(realpath ncbi-blast-2.12.0+/bin/update_blastdb.pl)



if [[ ! -z ${db_name} ]]; then
  echo "Downloading files for '${db_name}'"
  CMD="perl ${download_script} ${db_name}"
  echo ${CMD}
  eval ${CMD}
else
  echo "+-----------------------------+------------------------------------------------+"
  echo " File Name                    | Content Description                           "
  echo "+-----------------------------+------------------------------------------------+"
  echo "README                        | README for this subdirectory (this file)"
  echo "nr.*tar.gz                    | Non-redundant protein sequences from GenPept,"
  echo "                                Swissprot, PIR, PDF, PDB, and NCBI RefSeq"
  echo "nt.*tar.gz                    | Partially non-redundant nucleotide sequences from "
  echo "                                all traditional divisions of GenBank, EMBL, and DDBJ"
  echo "                                excluding GSS,STS, PAT, EST, HTG, and WGS."
  echo "landmark.tar.gz               | Proteome of 27 model organisms, see "
  echo "                                https://blast.ncbi.nlm.nih.gov/smartblast/smartBlast.cgi?CMD=Web&PAGE_TYPE=BlastDocs#searchSets"
  echo "16S_ribosomal_RNA             | 16S ribosomal RNA (Bacteria and Archaea type strains)"
  echo "18S_fungal_sequences.tar.gz   | 18S ribosomal RNA sequences (SSU) from Fungi type and reference material (BioProject PRJNA39195)"
  echo "28S_fungal_sequences.tar.gz   | 28S ribosomal RNA sequences (LSU) from Fungi type and reference material (BioProject PRJNA51803)"
  echo "ITS_RefSeq_Fungi.tar.gz       | Internal transcribed spacer region (ITS) from Fungi type and reference material (BioProject PRJNA177353)"
  echo "ITS_eukaryote_sequences.tar.gz| Internal transcribed spacer region (ITS) for eukaryotic sequences"
  echo "LSU_eukaryote_rRNA.tar.gz     | Large subunit ribosomal RNA sequences for eukaryotic sequences"
  echo "LSU_prokaryote_rRNA.tar.gz    | Large subunit ribosomal RNA sequences for prokaryotic sequences"
  echo "SSU_eukaryote_rRNA.tar.gz     | Small subunit ribosomal RNA sequences for eukaryotic sequences"
  echo "ref_euk_rep_genomes*tar.gz    | Refseq Representative Eukaryotic genomes (1000+ organisms)"
  echo "ref_prok_rep_genomes*tar.gz   | Refseq Representative Prokaryotic genomes (5700+ organisms)"
  echo "ref_viruses_rep_genomes*tar.gz   | Refseq Representative Virus genomes (9000+ organisms)"
  echo "ref_viroids_rep_genomes*tar.gz   | Refseq Representative Viroid genomes (46 organisms)"
  echo "refseq_protein.*tar.gz        | NCBI protein reference sequences"
  echo "refseq_rna.*tar.gz            | NCBI Transcript reference sequences"
  echo "swissprot.tar.gz              | Swiss-Prot sequence database (last major update)"
  echo "pataa.*tar.gz                 | Patent protein sequences"
  echo "patnt.*tar.gz                 | Patent nucleotide sequences. Both patent databases"
  echo "                                are directly from the USPTO, or from the EPO/JPO"
  echo "                                via EMBL/DDBJ"
  echo "pdbaa.*tar.gz                 | Sequences for the protein structure from the"
  echo "                                Protein Data Bank"
  echo "pdbnt.*tar.gz                 | Sequences for the nucleotide structure from the "
  echo "                                Protein Data Bank. They are NOT the protein coding"
  echo "                                sequences for the corresponding pdbaa entries."
  echo "taxdb.tar.gz                  | Additional taxonomy information for the databases "
  echo "                                listed here  providing common and scientific names"
  echo "FASTA/                        | Subdirectory for FASTA formatted sequences"
  echo "v4/                           | BLAST databases in version 4 (v4).  These files are no"
  echo "                                longer being updated."
  echo "cloud/	                      | Subdirectory of databases for BLAST AMI; see"
  echo "                                http://1.usa.gov/TJAnEt"
  echo "+-----------------------------+------------------------------------------------+"
  echo ""
  echo "For more info, see https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html"
  echo ""
  echo "Use this script to download databases"
  echo "perl ${download_script} <DB_NAME>"
  echo "  e.g."
  echo "    perl ${download_script} ref_euk_rep_genomes 	# Downloads all eukaryotic genomes"
fi