# DATA
Reference for where to retrieve sample data

## FASTQs
`get_fastqs.sh`
* Retrieves FASTQs from the 1000 Genomes Project. By default retrieves sample `ERR031940`, but can search for any files [1000 Genomes FTP Data Site](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/). Just go to the `/vol1/ftp/phase3/data/<SAMPLE>/sequence_read` path

```
$ ./get_fastqs.sh
RETRIEVING FQs
getting ERR031940_1.filt.fastq.gz
  done.
getting ERR031940_2.filt.fastq.gz
  done.
```

