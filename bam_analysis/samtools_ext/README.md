# Samtools Extension
Utils that seem like they could be in samtools

`extract_double_clipped_reads.sh` - Outputs `double_clipped_reads.sam` & `double_clipped_reads.fa` files with only the double-clipped reads from an input SAM

`isolate_clipped_reads.sh` - Outputs `*_clipped.sam` & `*_clipped_collated.bam` which are SAM/BAM files that have only the filtered double-soft/hard clipped reads.

## isolate_clipped_reads.sh
### Outputs
* `*_clipped.sam` - Reads w/ hard and/or soft clipping
* `*_clipped_collated.bam` - Collated `*_clipped.sam` where RGIDs are together
* `*_clipped_rgid.txt` - RGIDs matching input clipping filters

### Run
#### Extract all read pairs where both reads have double soft-clipping
```
$ ./isolate_clipped_reads.sh -f sample.sam -c S -p
[INPUTS]
BAM=sample.sam
CLIPPING_FILTER=S
PAIRED_CLIPPING_FILTER=TRUE

Extracting RGIDs of clipped BAM reads: S_clipped_rgid.txt
REGEX="[0-9]+[S].+[S]$"
Found 4329 clipped RGIDs

Including only RGIDs w/ clipping of both paired-reads
Found 2003 RGIDs with clipping of both pairs

Writing sample_clipped.sam
	headers...
	reads...

Shuffling/Grouping reads into collated BAM=sample_clipped_collated.bam

Done.
```
#### Extract all reads pairs where at least one of the reads have double soft-clipping
```
$ ./isolate_clipped_reads.sh -f sample.sam -c H
[INPUTS]
BAM=sample.sam
CLIPPING_FILTER=H
PAIRED_CLIPPING_FILTER=

Extracting RGIDs of clipped BAM reads: H_clipped_rgid.txt
REGEX="[0-9]+[H].+[H]$"
Found 3 clipped RGIDs

Including all RGIDs w/ clipping of at least one paired-read

Writing sample_clipped.sam
	headers...
	reads...

Shuffling/Grouping reads into collated BAM=sample_clipped_collated.bam

Done.
```

#### Extract all reads pairs where both reads have double soft OR hard clipping
```
$ ./isolate_clipped_reads.sh -f sample.sam -c HS -p
[INPUTS]
BAM=sample.sam
CLIPPING_FILTER=H
PAIRED_CLIPPING_FILTER=

Extracting RGIDs of clipped BAM reads: H_clipped_rgid.txt
REGEX="[0-9]+[H].+[H]$"
Found 3 clipped RGIDs

Including all RGIDs w/ clipping of at least one paired-read

Writing sample_clipped.sam
	headers...
	reads...

Shuffling/Grouping reads into collated BAM=sample_clipped_collated.bam

Done.
[pplnuser@sga01 test_3]$ ../isolate_clipped_reads.sh -f sample.sam -c HS -p
[INPUTS]
BAM=sample.sam
CLIPPING_FILTER=HS
PAIRED_CLIPPING_FILTER=TRUE

Extracting RGIDs of clipped BAM reads: HS_clipped_rgid.txt
REGEX="[0-9]+[HS].+[HS]$"
Found 4332 clipped RGIDs

Including only RGIDs w/ clipping of both paired-reads
Found 2003 RGIDs with clipping of both pairs

Writing sample_clipped.sam
	headers...
	reads...

Shuffling/Grouping reads into collated BAM=sample_clipped_collated.bam

Done.
```

