# Blast Server Setup

## Description
Automates downloading of blast DB databases from scratch. The process is already made very simple by [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html), but `setup_blast.sh` is meant to set everything up in a single command.

### NCBI
NCBI already has pre-formatted databases and a convenient downloadable script, `update_blastdb.pl`, available via FTP. In general, the steps are - 
  1. Download and extract all preformatted databases with the same prefix from NCBI's FTP server (or make own) - *Note - this can take hours*
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
# ... come back in a few hours
```

### Inputs
* `-v`: `string` (optional), blast+ version to download. See ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+ for more details, *Default: Latest available*
  > `2.2.18`-`2.12.0`, `l` (for latest)
* `-o`: `string` (optional) - OS version for blast+ scripts, *Default: `x64-linux`*
  > `win64`, `x64-linux`, `x64-macosx`, `x64-win64`
* `-d`: `string` (optional) -  `update_blastdb.pl` downloads all `*tar.gz` w/ this prefix in ftp://ftp.ncbi.nlm.nih.gov/blast/db. For a full list, see [blastdb README](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html)
  > **SINGLE**: `ref_euk_rep_genomes`, `nt`, `nr`, `ref_prok_rep_genomes`, `ref_viruses_rep_genomes`, `ref_viroids_rep_genomes`, `refseq_protein`, `refseq_rna`

  > **MULTIPLE**: `-d "human_genome ref_prok_rep_genomes nt ref_euk_rep_genomes"`
* `-p`: `string` (optional, default: current working directory) - If specifying a database to download, `-p` specifies the path to write those files to. Note - this will later be the directory pointed to by the environment variable, `BLASTDB`, when running blast. See [Configuring BLAST](https://www.ncbi.nlm.nih.gov/books/NBK569858/)
  > `-p ~/blast_db/preformatted_dbs`


## Manual Download Guide
### FTP Notes for Manual Download
These are paraphrased from the [blastdb README](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html).
* **Pre-formatted Databases** Use the pre-formatted databases, or follow the README to create your own
  * "Pre-formatted databases must be downloaded using the update_blastdb.pl script or via FTP in binary mode" ([REF](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html))
* **Extraction** If using the pre-formatted databases, they must be untarred (`tar -zxvf *.tar.gz`) before use.
  * As a note, I found that the extracted DBs are not much larger than their tar'd versions. For instance, when I downloaded `ref_euk_rep_genomes`, the tar'd was ~240GB and the untar'd was ~250GB.
  * **Issue with `--decompress` option** The `--decompress` option for `update_blastdb.pl` is intended to "decompresses the archives in the current working directory...and delete(s) the downloaded archive to save disk space" ([REF](https://www.ncbi.nlm.nih.gov/books/NBK62345/)). However I always receive the error below *with at least 35 GB of disk space left* so I removed this option and manually added extraction and deletion of each file in. `setup_blast.sh`. It seems this is what the `--decompress` option would have done anyway because **all the tar files are downlaoded** before any extraction lines are logged. I assume each `*.tar.gz` file will be extracted and deleted individually (rather than save all deletions until the end), which is what `setup_blast.sh` does.

      *Error for `--decompress` option*
      ```
      Could not write data to './preformatted_dbs/ref_euk_rep_genomes.00.nsq' at ./ncbi-blast-2.12.0+/bin/update_blastdb.pl line 496.
      Failed to decompress ref_euk_rep_genomes.00.tar.gz (Could not write data to './preformatted_dbs/ref_euk_rep_genomes.00.nsq'), please do so manually.
      ```

### Perl Script
1. Determine desired blast version & os (e.g. "2.12.0" and "linux-64" respectively)
2. Download NCBI's tar file from the appropriate FTP folder at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast%2B (e.g. `curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.12.0/ncbi-blast-2.12.0+-x64-linux.tar.gz -o ncbi-blast-2.12.0+-x64-linux.tar.gz`)
    ```
    version=2.12.0
    os=linux-64
    curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${version}/ncbi-blast-${version}+-${os}.tar.gz
    ```
  * NOTE - Large databases, which "[are formatted in multiple one-gigabyte volumes](https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html)",  can take hours. E.g. ref_euk_rep_genomes is 90 files, each about 3GB, and can take about 3-4 hours at 20 MBPS download 
   
3. Extract tar
4. Run the `update_blastdb.pl` w/ the prefix of the FASTA files you want
    ```
    perl ncbi-blast-2.12.0+/bin/update_blastdb.pl ref_euk_rep_genomes
    ```
    
### curl/wget
* For some reason, sometimes the perl script doesn't download the `.md5` files and the entire script fails
```
DB=nt # One of the preformatted names
md5sum_calc="md5sum_calc.txt"
md5sum_download="md5sum_dwld.txt"
db_files=$(curl ftp://ftp.ncbi.nlm.nih.gov/blast/db/ 2> /dev/null | grep -oE "\s${DB}.[0-9]+.tar.gz" | cut -d' ' -f2)
for file in ${db_files}; do
  curl ftp://ftp.ncbi.nlm.nih.gov/blast/db/${file} > ${file}
  if [[ $? -eq 0 ]]; then
    echo "SUCCESS=${file}"
  else
    echo "FAIL=${file}"
  fi
  md5_file=${file}.md5
  curl ftp://ftp.ncbi.nlm.nih.gov/blast/db/${md5_file} > ${md5_file}
  md5sum ${file} >> ${md5sum_calc}
  cat ${md5_file} >> ${md5sum_download}
done

diff ${md5sum_calc} ${md5sum_download} # Should Match & not have any diffs
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
**STEP 0** - `blastdbcheck`
* In the ncbi executable download, there is a script that checks the preformatted databses (`./ncbi-blast-2.12.0+/bin/blastdbcheck`). Run this and verify there are no errors

#### Notes on Specific Errors
* `Error: mdb_env_open` - Re-download and extract preformatted databases
* `BLAST Database error: Cannot memory map` - Not sure, but most likely a resource issue. If there is a line, stating `Number of files opened: ###` and that number is less than the total files in the `BLASTDB` directory for the DB you're using, it might be that your system has a limit on the number of open files. See below,

  ```
  $ ulimit -f -n
  file size               (blocks, -f) unlimited
  open files                      (-n) 256
  ```
  Check if "open files" is less than the total files in the untarred database, if it is, this needs to be changed. Look up "modify limit of file descriptors" for your OS and change it to a 2^x greater than the number of files (e.g. 8192). e.g. For macOS,
  ```
  $ sudo launchctl limit maxfiles 8192 unlimited  # CHANGE
  $ ulimit -n 8192
  $ launchctl limit maxfiles                      # VERIFY
  	maxfiles    8192           10240
  $ ulimit -f -n
    file size               (blocks, -f) unlimited
    open files                      (-n) 8192
  # NOTE - this is a temporary change until the MAC is rebooted or you log out
  ```
* Nothing is printed in BLAST results
  * Probably a memory issue
  * Try to see if some results are printed by temporarily changing the environment variable `BATCH_SIZE` to something small. This is the number of queries that will be concatenated to conserve memory. E.g.
  ```
  export BATCH_SIZE=10
  ```
  * As a note, there are also additional environment variables used by blast. See [Configuring BLAST via environment variables](https://www.ncbi.nlm.nih.gov/books/NBK569858/#_usrman_Config_BLAST_Configuring_BLAST_vi_1_)
* `Critical: Failed to initialize SSL provider MBEDTLS: Unknown` - Not sure, but maybe related to fire wall. See below,
  * [NCBI Firewall Info](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/NETWORK/firewall.html)
  * [Check firewall ports](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/NETWORK/fwd_check.cgi) - Sometimes blast needs to query NCBI even when running locally, e.g. when running w/ `-remote`
* **No error, and no results** - blast is a cpu/memory-intensive process. It is not uncommon to wait hours and even days for large queries. Try these debug steps - 

  1. Check that it is still running
    ```
    $ ps -ef | grep blast   # Should get more than the `grep blast` command
    502 41870 41867   0 11:03AM ttys002   13:21.40 ./ncbi-blast-2.12.0+/bin/blastn -db ref_euk_rep_genomes -num_threads 7 -outfmt 6 -query query.fa -out results.out
    502 42597 35172   0 11:10AM ttys002    0:00.00 grep blast
    ```
  2. If multiple queries are being run in the query file, try one and increase the number of threads, `-num_threads` option. BLAST will [query pack](https://doctorlib.info/medical/blast/13.html) with a set of queries and only return results when ALL of them have been analyzed. 
  3. Re-run w/ an optimized query

#### Other
* There seems to be a character limit for `blastn -db ${DBS}...` - **TODO**

### References
* [blastn documentation](https://www.ncbi.nlm.nih.gov/books/NBK569856/)

