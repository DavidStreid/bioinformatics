#!/bin/bash
tax_id=$1
name=$2

DB=ref_viruses_rep_genomes

fname="$(echo ${name} | sed 's/ /_/g')__${tax_id}.fa"
blastdbcmd -db ${DB} -taxids ${tax_id} > ${fname}
echo ${fname}