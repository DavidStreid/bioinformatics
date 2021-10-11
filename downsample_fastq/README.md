# Downsample FASTQ file

## Setup
```
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make
# Add seqtk to PATH
```
See [seqtk documentation](https://github.com/lh3/seqtk#introduction)

## Run
Usage: `./downsample.sh -f [input_bam] -r [num_reads]`

```
# Downsample FASTQ to 1000 reads
./downsample.sh -f sample.fastq -r 1000
```
