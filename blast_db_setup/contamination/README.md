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
$ python3 analyze.py -f blast_results.tsv -p my_sample
SETTINGS
	MAX_E_VALUE=0.1
	MAX_MAGNITUDE_DIFFERENCE_ALLOWED=1
[IN] BLAST RESULT SUMMARY
	INPUT = blast_results.tsv
	Total Reads = 93
AGGREGATION
	total_blast_results=843597
	Total Reads with Valid Blast Results = 93
	MOST IDENTIFIED=firmicutes (0.68)
[OUT] RGID to Identification:	my_sample___query_identities.tsv
[OUT] ID Proportion Summary:	my_sample___identity_summary.tsv
[OUT] Pie Chart (Total IDs=4):	my_sample___blast_identifications_per_read_id.pdf
[OUT] BLAST RESULT SUMMARY:	my_sample___blast_result_analysis.tsv
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