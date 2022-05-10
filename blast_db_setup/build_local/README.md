# Local BLASTDB Setup
* This provides a script, `taxid_to_blastdb.sh`, to go from taxonomic ID to ready-to-go BLASTDB files created from NCBI's `makeblastdb` script

## Prerequisites
* Have the `makeblastdb` executable in PATH. To do this, download a blast package via [`setup_blast.sh`](https://github.com/DavidStreid/bioinformatics/blob/main/blast_db_setup/setup_blast.sh).

## Run
Two scripts will create BLAST DBs from input taxonomic IDs. The difference is in where the reference files are taken from.
### Extract all fastas for the given taxonomic ID from the input BLAST DB
```
$ ./extract_taxid_blastdb.sh 8036
TAXID=8036
SPECIES=Salvelinus_alpinus
BLAST_DB=nt
8036___Salvelinus_alpinus___nt
        preparing taxid fasta: 8036___Salvelinus_alpinus.fa
\tRemoving spaces from fasta description lines
        Preparing tax ID map: taxid_map__Salvelinus_alpinus__8036.txt
        Creating BLAST DB: Salvelinus_alpinus__8036
        makeblastdb -in /Users/dstreid/work/repos/bioinformatics/blast_db_setup/build_local/8036___Salvelinus_alpinus___nt/fasta/8036___Salvelinus_alpinus.fa -parse_seqids -taxid_map /Users/dstreid/work/repos/bioinformatics/blast_db_setup/build_local/8036___Salvelinus_alpinus___nt/fasta/taxid_map__Salvelinus_alpinus__8036.txt -title 'Salvelinus_alpinus__8036' -dbtype nucl -out Salvelinus_alpinus__8036
        SUCCESS TAXID=8036 SPECIES=Salvelinus_alpinus
```

### Take from NCBI's mapping of taxnomic ID
The `taxid_to_blastdb.sh` script will take in a list of taxonomic id arguments. It will then create a directory contatining the source fastas and BLASTDB files if a valid taxonomic ID to fasta file exists for that taxonomic ID.
```
$ TAX_ID_LIST="889201 1897061"
$ ./taxid_to_blastdb.sh ${TAX_ID_LIST}
Locating sources for list of taxonomic IDs
Finding source for 889201
...
889201  Streptococcus cristatus ATCC 51100      na      Scaffold        https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/187/855/GCA_000187855.1_ASM18785v1/GCA_000187855.1_ASM18785v1_genomic.fna.gz
Finding source for 1897061
...
1897061 Cryobacterium sp. SO1   na      Scaffold        https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/210/215/GCA_004210215.1_ASM421021v1/GCA_004210215.1_ASM421021v1_genomic.fna.gz
Creating BLASTDB files
1897061___Cryobacterium_sp_SO1___na___Scaffold
        downloading https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/210/215/GCA_004210215.1_ASM421021v1/GCA_004210215.1_ASM421021v1_genomic.fna.gz
        ...
        SUCCESS TAXID=1897061 SPECIES=Cryobacterium_sp_SO1
Done.
```

## Background
* The preformatted NCBI BLAST databases are very comprehensive, but can be too large to store in memory. Selecting a few organisms to build a database from can make running BLAST more efficient and effective. This is done via the `makeblastdb` script and an input reference file
* Inputs:
  * FASTA file (required)
  * TaxIDs, i.e "Taxonomic node": An ID that identifies sequences in the FASTA file (these can be returend in the results when running BLAST
    * These are basically the description lines of the fasta file
    * No spaces or `|` characters

### Run `makeblastdb`
* See [Building a BLAST database with your (local) sequences](https://www.ncbi.nlm.nih.gov/books/NBK569841/) for more information
```
$ head -3 ${tax_id_map}   # input fasta sequence IDs to their taxonomic node (not mandatory, but highly recommended)
seq1 9501
seq2 812
seq3 10001
$ makeblastdb -in ref_genome.fa -parse_seqids -taxid_map ${tax_id_map} -dbtype nucl      # Basic blastn local BLAST DB setup
```
  
### Find a Reference
#### From taxonomic ID
NCBI provides a mapping of taxonomic IDs to their reference/representative genomes, or any assembly related to that ID. This file is available [here](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt).
The `staxis_to_ncbi_assemlby.sh` is provided to take an NCBI ID and provide a link to the best assembly to use as input to the `makeblastdb` script.
##### Run
This will run just the portion that will identify a valid reference file to use given a taxonomic ID
```
$ ./staxis_to_ncbi_assemlby.sh 9606 2> err.out    # Returns the reference genome for Homo Sapiens (Taxonomy ID=9606)
9606    reference genome        Chromosome      https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz
```
