# Blast Server Setup

## Description
Automates downloading of blast DB databases from scratch. The process is already made very simple by [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html), but `setup_blast.sh` is meant to set everything up in a single command.

### NCBI
NCBI already has pre-formatted databases and a convenient downloadable script, `update_blastdb.pl`, available via FTP. But, to avoid navigating the FTP repo ,`ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+`, and all its versions for this (which is just downloading the right `*tar.gz` for the OS and picking a version), I plan to use `setup_blast.sh` in the future to make my life easier.

If not using the script, here are some notes paraphrased from the `blastdb` [README](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html).
* Use the pre-formatted databases, or follow the README to create your own
* "Pre-formatted databases must be downloaded using the update_blastdb.pl script or via FTP in binary mode" [REF](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html)
* If using the pre-formatted databases, they must be untarred (`tar -zxvf *.tar.gz`) before use. As a note, I found that the extracted DBs are not much larger than their tar'd versions. For instance, when I downloaded `ref_euk_rep_genomes`, the tar'd was ~240GB and the untar'd was ~250GB.

## Run
### Inputs
* `-v`: `string` (optional), blast+ version to download. See ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+ for more details
  > `2.2.18` to `2.12.0`
* `-o`: `string` (optional) - OS version for blast+ scripts
  > `win64`, `x64-linux`, `x64-macosx`, `x64-win64`
* `-d`: `string` (optional) -  `update_blastdb.pl` downloads all `*tar.gz` w/ this prefix in ftp://ftp.ncbi.nlm.nih.gov/blast/db. For a full list, see [blastdb README](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html)
  > `ref_euk_rep_genomes`, `nt`, `nr`
* `-p`: `string` (optional, default: current working directory) - If specifying a database to download, `-p` specifies the path to write those files to. Note - this will later be the directory pointed to by the environment variable, `BLASTDB`, when running blast. See [Configuring BLAST](https://www.ncbi.nlm.nih.gov/books/NBK569858/)
  > `-p ~/blast_db/preformatted_dbs`

### Download script, then download DB later
```
./setup_blast.sh
```

### Download everything including DB 
```
version=2.12.0
os=linux-64
db_name=ref_euk_rep_genomes

./setup_blast.sh -v ${version} -o ${x64-linux} -d ${db_name}
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

## Run blast
```
# Assign this environment variable to the directory w/ the db files that have been extracted,
#  i.e. directory w/ `*.nhr`, `*.nin`, `*.nnd`, `*.nnd`, `*.nni`, `*.nog`, `*.nsq`
BLASTDB=...       
DB_NAME=ref_euk_rep_genomes          # Name of the database to use (prefix of the DB file names)
./ncbi-blast-2.12.0+/bin/blastn -db ${DB_NAME} -query sample.fa -out results.out
```

### Troubleshooting
* `Error: mdb_env_open` - Re-download and extract preformatted databases
