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
    * The chosen database MUST already downloaded and available in the `BLASTDB` directory. The table below contains databases I think are good choices for detecting contamination

		|Name|Type|Description|
		|:---:|:---:|:---:|
		|nt|DNA|Comprensive set of sequences from GenBank, EMBL, DDBJ, PDB, & RefSeq (If only choosing one - choose this one)|
		|ref_euk_rep_genomes|DNA|Representative Eukaryotic Genomes from NCBI's Refseq Genomes database|
		|ref_prok_rep_genomes|DNA|Representative Eukaryotic Genomes from NCBI's Refseq Genomes database|
		|ref_viruses_rep_genomes|DNA|Representative Virus (not viroid) Genomes from NCBI's Refseq Genomes database|
		|human_genome|DNA|Homo sapiens GRCh38|
		|env_nt|DNA|samples taken directly from their environment (e.g. Mine Drainage projects)|

    * Not every read is guaranteed to return results based on the default criteria and database used. Therefore, choosing the `-d` database parameter is crucial for getting accurate results
	    * For instance, If the DB used does not include the species, but is large enough (note: E-Value is directly related to database size) that the read maps to a highly similar region of another species in the DB, then a BLAST result could return a seemingly good E-Value, say `E-Value=1e-10`, that is actually incorrect. In this hypothetical, adding a database with the correct species could reduce the E-Value several orders of magnitude to something like `E-Value=1e-50`, which would clearly be the more confident identification.



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

#### References
* [Github Docs](https://github.com/ncbi/blast_plus_docs#blast-databases)
* [BLAST DB Help](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html)

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