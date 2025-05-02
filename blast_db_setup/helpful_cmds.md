## Retrieve Sequences
```
$ blastdbcmd -db ref_viruses_rep_genomes -taxids 11137
>NC_002645.1 Human coronavirus 229E, complete genome
ACTTAAGTACCTTATCTATCTACAGATAGAAAAGTTGCTTTTTAGACTTTGT...
```

## Available Databases for Download

* `DB` used later

```
DB=ref_euk_rep_genomes # ref_euk_rep_genomes nt nr ref_prok_rep_genomes ref_viruses_rep_genomes ref_viroids_rep_genomes refseq_protein refseq_rna
```

## Nucleotide BLAST
```
blastn -db ${DB} -query query.fa -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue"
```

## Create Database
* Taxonomic ID doesn't seem to like spaces in the name, or special characters

```
FASTA=$1
OUT=$2        # path_to/db_name_prefix
TAXID_MAP=$3

# $ cat ${TAXID_MAP}
# seq1 9606 # sequence 1 maps to taxonomic ID for human (9606)

CMD="makeblastdb -in ${FASTA} -taxid_map tax_id.txt -parse_seqids -dbtype nucl -out ${OUT}"
```

