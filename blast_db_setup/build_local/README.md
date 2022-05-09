## Prerequisites
* Download a blast package via `setup_blast.sh`

## Info
* The preformatted NCBI BLAST databases are very comprehensive, but can be too large to store in memory. Selecting a few organisms to build a database from can make running BLAST more efficient and effective. This is done via the `makeblastdb` script and an input reference file
* References
  * [Building a BLAST database with your (local) sequences](https://www.ncbi.nlm.nih.gov/books/NBK569841/)
  
## Run
```
$ makeblastdb -in ref_genome.fa -dbtype nucl      # Basic blastn local BLAST DB setup
```

## Find a Reference
### From taxonomic ID
NCBI provides a mapping of taxonomic IDs to their reference/representative genomes, or any assembly related to that ID. This file is available [here](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt).
The `staxis_to_ncbi_assemlby.sh` is provided to take an NCBI ID and provide a link to the best assembly to use as input to the `makeblastdb` script.
#### Run
```
./staxis_to_ncbi_assemlby.sh 9606     # Returns the reference genome for Homo Sapiens (Taxonomy ID=9606)
```
