# Blast Result analyzer
Ouputs a processed TSV of the best result for each read along w/ helpful statistics such as `magnitude_of_best_result` & `proportion_valid_blast_results`, which allow for analysis of what results are valid.
* e.g. If `magnitude_of_best_results` is very high, this is probably a valid hit

## Run
### Create BLAST results
```
$ ./blast_fasta.sh sample.fa
```
* Note - This is set to blast against the eukaryotic `ref_euk_rep_genomes`. To modify this, change `blast_fasta.sh` like below to one of the databases available on the [NCBI website](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html) (assuming this database is already downloaded and is available in the `BLASTDB` location)
```
...
DB=ref_prok_rep_genomes 
...
```


### Find best identification for each read that returns results
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

## Outputs
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