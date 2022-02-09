# Blast Server Setup

## Description
Automates downloading of blast DB databases from scratch. The process is already made very simple by [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html), but `setup_blast.sh` is meant to set everything up in a single command.

### NCBI
NCBI already has pre-formatted databases and a convenient downloadable script, `update_blastdb.pl`, available via FTP. In general, the steps are - 
  1. Download and extract all preformatted databases with the same prefix from NCBI's FTP server (or make own) - *Note - this usually takes hours*
  2. Download a compatible executable of the desired blast from NCBI 
  3. Run blast specifying the prefix of the downloaded files of step 1 as the database to blast against

See [Manual Download Guide](#manual-download-guide) for more details. But, to avoid navigating the FTP repo ,`ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+`, and all its versions, there's `setup_blast.sh` that does steps 1 & 2.

## Run
### Download executables with defaults
```
./setup_blast.sh
```

### Download executables specifying version
```
version=2.11.0
os=x64-macosx

./setup_blast.sh -v ${version} -o ${os}
```

### Download executables and preformatted db files with defaults
```
db_name=ref_euk_rep_genomes
blastdb=./preformatted_db

./setup_blast.sh -d ${db_name} -p ${blastdb}      # Could specify different version & os w/ -v & -o 
```

### Inputs
* `-v`: `string` (optional), blast+ version to download. See ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+ for more details, *Default: Latest available*
  > `2.2.18`-`2.12.0`, `l` (for latest)
* `-o`: `string` (optional) - OS version for blast+ scripts, *Default: `x64-linux`*
  > `win64`, `x64-linux`, `x64-macosx`, `x64-win64`
* `-d`: `string` (optional) -  `update_blastdb.pl` downloads all `*tar.gz` w/ this prefix in ftp://ftp.ncbi.nlm.nih.gov/blast/db. For a full list, see [blastdb README](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html)
  > `ref_euk_rep_genomes`, `nt`, `nr`, `ref_prok_rep_genomes`, `ref_viruses_rep_genomes`, `ref_viroids_rep_genomes`, `refseq_protein`, `refseq_rna`
* `-p`: `string` (optional, default: current working directory) - If specifying a database to download, `-p` specifies the path to write those files to. Note - this will later be the directory pointed to by the environment variable, `BLASTDB`, when running blast. See [Configuring BLAST](https://www.ncbi.nlm.nih.gov/books/NBK569858/)
  > `-p ~/blast_db/preformatted_dbs`


## Manual Download Guide
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
    
### FTP Notes for Manual Download
These are paraphrased from the [blastdb README](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html).
* **Pre-formatted Databases** Use the pre-formatted databases, or follow the README to create your own
  * "Pre-formatted databases must be downloaded using the update_blastdb.pl script or via FTP in binary mode" ([REF])(https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html)
* **Extraction** If using the pre-formatted databases, they must be untarred (`tar -zxvf *.tar.gz`) before use.
  * As a note, I found that the extracted DBs are not much larger than their tar'd versions. For instance, when I downloaded `ref_euk_rep_genomes`, the tar'd was ~240GB and the untar'd was ~250GB.
  * **Issue with `--decompress` option** The `--decompress` option for `update_blastdb.pl` is intended to "decompresses the archives in the current working directory...and delete(s) the downloaded archive to save disk space" ([REF](https://www.ncbi.nlm.nih.gov/books/NBK62345/)). However I always receive the error below *with at least 35 GB of disk space left* so I removed this option and manually added extraction and deletion of each file in. `setup_blast.sh`. It seems this is what the `--decompress` option would have done anyway because **all the tar files are downlaoded** before any extraction lines are logged. I assume each `*.tar.gz` file will be extracted and deleted individually (rather than save all deletions until the end), which is what `setup_blast.sh` does.

      *Error for `--decompress` option*
      ```
      Could not write data to './preformatted_dbs/ref_euk_rep_genomes.00.nsq' at ./ncbi-blast-2.12.0+/bin/update_blastdb.pl line 496.
      Failed to decompress ref_euk_rep_genomes.00.tar.gz (Could not write data to './preformatted_dbs/ref_euk_rep_genomes.00.nsq'), please do so manually.
      ```


## Run blast
Once the `bin` of all the bash executables for whatever version was specified and the pre-formatted database files are downloaded, blast is very simple to run,
1. Export `BLASTDB` to be the directory w/ the downloaded database files. See [Configuring BLAST](https://www.ncbi.nlm.nih.gov/books/NBK569858/) for more information. This is the directory w/ the db files that have been extracted, i.e. directory w/ `*.nhr`, `*.nin`, `*.nnd`, `*.nnd`, `*.nni`, `*.nog`, `*.nsq`
    ```
    $ export BLASTDB=...
    ```
2. Run the blast executable of choice, specifying the desired database to try to map to
    ```
    $ ./ncbi-blast-2.12.0+/bin/blastn -db ref_euk_rep_genomes -query sample.fa -out results.out
    ```

**Example**
```
export BLASTDB==path/to/formatted_db_files
DB_NAME=ref_euk_rep_genomes          # Name of the database to use (prefix of the DB file names)
./ncbi-blast-2.12.0+/bin/blastn -db ${DB_NAME} -query sample.fa -out results.out
```

### Troubleshooting
* `Error: mdb_env_open` - Re-download and extract preformatted databases
* `BLAST Database error: Cannot memory map` - Memory issue, need to run on a server w/ more memory
### References
* [blastn documentation](https://www.ncbi.nlm.nih.gov/books/NBK569856/)

