# Helpers for blastn analysis
## Run
### Run BLAST
BLASTs an input fasta against a specified database
```
$ ./blast_fasta.sh -f sample.fa [-d blastn database]
```
* `f`: string (required), fasta file
> Eg: `sample.fa`
* `d`: string, previously downloaded blastn database to use See [BLAST DB Help](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html) for list of databases
> Eg: `ref_prok_rep_genomes`, `nt`, `"ref_prok_rep_genomes nt"`
  * Default: `nt` (this is the [blastn default](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) and combines nucleotide sequences from across GenBank, EMBL, and DDBJ)
  * Notes:
    * This database MUST already downloaded and available in the `BLASTDB` directory
    * Not every read is guaranteed to return results based on the default criteria and database used

#### Outputs
`blast_results.tsv` - TSV w/ the following fields for each read that returns BLAST results
* `qaccver` - RGID
* `saccver `
* `pident`
* `length` 
* `mismatch` 
* `gapopen` 
* `qstart` 
* `qend` 
* `sstart` 
* `send` 
* `evalue` - likelihood that the BLAST result would have been returned by chance
* `bitscore` 
* `staxids`
* `sscinames` 
* `scomnames` - Common name for BLAST result 
* `sblastnames` 
* `sskingdoms` - Kingdom of BLAST result

### Extract best identifications for each read
Ouputs a processed TSV of the best result for each read along w/ helpful statistics such as `magnitude_of_best_result` & `proportion_valid_blast_results`, which allow for analysis of what results are valid.
* e.g. If `magnitude_of_best_results` is very high, this is probably a valid hit
```
$ python3 analyze.py -f blast_results.tsv
SETTINGS
        MAX_E_VALUE=0.1
        MAX_MAGNITUDE_DIFFERENCE_ALLOWED=3
[IN] BLAST RESULT SUMMARY
        INPUT = blast_results.tsv
        Total Reads = 1095
AGGREGATION
        total_blast_results=784379
        Total Reads with Valid Blast Results = 1095
IDENTIFICATIONS
        primates: 1052
        even-toed ungulates: 14
        carnivores: 2
        placentals: 1
        eudicots: 1
        ants: 5
        apicomplexans: 11
        odd-toed ungulates: 1
        eukaryotes: 3
        bats: 2
        monocots: 1
        coelacanths: 1
        marsupials: 1
MOST_IDENTIFIED=primates
IDENTIFICATION_PROPORTION=0.960730593607306
[OUT] GRAPHING
        file=blast_results.pdf (n=13)
[OUT] BLAST RESULT SUMMARY
        IDS & STATISTICS=output.tsv
        IDS ONLY=query_identities.tsv
```

#### Outputs
`blast_results.pdf` - Pie Chart of Identifies

`query_identities.tsv` - TSV of just identifies w/o statistics

`identity_summary.tsv` - Identifications and their proportions

`blast_result_analysis.tsv` - TSV of identifies & statistics for analysis
* qaccver
* pident
* evalue
* magnitude_of_best_result
* proportion_valid_blast_results
* bitscore
* scomnames
* sblastnames
* sskingdoms