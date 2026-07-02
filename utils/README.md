# Utilities

## Description

General-purpose tools that come up in bioinformatic processing

## Tools

### `gzip -t` - checks the compressed gzip'd file's integrity

* no output and exitcode=0 if no issues
* output and non-zero exitcode

Checks whether a gzip file is valid, e.g. FASTQ (`*.fastq.gz`)

```
$ gzip -t valid.fastq.gz
$ gzip -t invalid.fastq.gz
gzip: invalid.fastq.gz: unexpected end of file
```
