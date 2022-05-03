# Samtools Extension
Utils that seem like they could be in samtools

`extract_double_clipped_reads.sh` - Outputs `double_clipped_reads.sam` & `double_clipped_reads.fa` files with only the double-clipped reads from an input SAM
`isolate_clipped_reads.sh` - Outputs `*_clipped.sam` & `*_clipped_collated.bam` which are SAM/BAM files that have only the filtered double-soft/hard clipped reads.
* `*_clipped.sam` - Only reads with clippings
* `*_clipped_collated.bam` - All reads where one of the read pairs has clipping


## isolate_clipped_reads.sh
```
./isolate_clipped_reads.sh -b sample.bam -f S
BAM=sample.bam
FILTER=S
	filter - H=Hard	S=Soft	HS=Hard & Soft
Extracting RGIDs of clipped BAM reads: soft_clipped_rgid.txt
REGEX=[0-9]+[S].+[S]$
Found 12 clipped RGIDs
Writing partial_clipped.sam
	headers...
	reads...
Writing Collated BAM=sample_clipped_collated.bam
Done.
```