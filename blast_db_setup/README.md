# Blast Server Setup

## Description
Notes for downloading/updating a BLAST database locally. Most important thing is getting the `update_blastdb.pl`

Inlcudes script that (hopefully) automates this process (`setup_blast.sh`)

## Run
### Inputs
* `version`: `string` (optional), blast+ version to download. See ftp.ncbi.nlm.nih.gov/blast/executables/blast+
  > `2.2.18` to `2.12.0`
* `os`: `string` (optional) - OS version for blast+ scripts
  > `win64`, `x64-linux`, `x64-macosx`, `x64-win64`
* `dbname`: `string` (optional) - `update_blastdb.pl` downloads all `*tar.gz` w/ this prefix in ftp://ftp.ncbi.nlm.nih.gov/blast/db. For a full list, see [blastdb README](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html)
  > `ref_euk_rep_genomes`, `nt`, `nr`

### Download script, then download DB later
```
./setup_blast.sh
```

### Download everything including DB 
```
version=2.12.0
os=linux-64
db_name=ref_euk_rep_genomes

./setup_blast.sh ${version} ${os} ${db_name}
```



## Run (w/o `setup_blast.sh`)
1. Determine desired blast version & os (e.g. "2.12.0" and "linux-64" respectively)
2. Download NCBI's tar file from the appropriate FTP folder at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast%2B (e.g. `curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.12.0/ncbi-blast-2.12.0+-x64-linux.tar.gz -o ncbi-blast-2.12.0+-x64-linux.tar.gz`)
```
version=2.12.0
os=linux-64
curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${version}/ncbi-blast-${version}+-${os}.tar.gz
```
3. Extract tar
4. Run the `update_blastdb.pl` w/ the prefix of the FASTA files you want
```
perl ncbi-blast-2.12.0+/bin/update_blastdb.pl ref_euk_rep_genomes
```
