# Simple
Simple things 

## Estimate alignment time - BWA MEM will output the number of processed reads 

1. Get Rate, **R** - reads processed / sec - e.g. `67742 reads / 41.871 sec = 1618 reads/sec.`
```
$ grep -A 1 -B 1 -E "Processed [0-9]+ reads"
[0000] 10. Calling kt_for - worker_sam
	[0000][ M::mem_process_seqs] Processed 67742 reads in 42.046 CPU sec, 41.871 real sec
[0000] 2. Calling mem_process_seqs.., task: 5046
```

2. Get Total Reads, **T**
```
$ zcat ${FQ} | wc -l
# OR
$ seqkit stats ${FQ}
```

3. Time Estimate, **E** - `E = R / T`

## Sort, compress, & index VCF
1. Sort
```
$ java -jar picard.jar SortVcf \
      I=unsorted.vcf \
      O=sorted.vcf
```

2. Compress - [`bgzip`](http://www.htslib.org/doc/bgzip.html) provides a way that allows indexes to be built on the compressed file
```
$ bgzip sorted.vcf
```

3. Index - [`tabix`](http://www.htslib.org/doc/tabix.html)
```
$ tabix -p vcf sorted.vcf.gz
```

4. Review
```
$ ls
sorted.vcf.gz
sorted.vcf.gz.tbi
```

## Get Total Number of bases in BED file
```
cat file.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
```
* REF - [BioStars](https://www.biostars.org/p/68283/#68292)

## Extracting R1/R2 FASTQs from BAMs
* Rembmer to sort BAM by name so that paired reads aren't missed

```
echo "Sorting SAM=${sorted_sam}"
samtools sort -n -o ${sorted_sam} ${original_sam}
echo "Extracting reads in FQ=${fq1} FQ2=${fq2}"
bedtools bamtofastq -i ${sorted_sam} -fq ${fq1} -fq2 ${fq2}
```

### Verify
```
samtools flagstat ${original_sam}
seqkit stats ${fq1}
seqkit stats ${fq2}
```
*  INPUT == OUTPUT
  * If paired,
    * INPUT (seqkit) = num_seqs (R1) + num_seqs (R2)
    * OUTPUT (samtools) = paired in sequencing
  * If single-read,
    * INPUT (seqkit) = num_seqs (R1)
    * OUTPUT (samtools) = mapped
  *  [REF](https://github.com/arq5x/bedtools2/issues/797)

## SAM/BAM

### samtools -F/-f
**Description** - Use SAM flags to filter on/out reads
* **`-F`** - Don't Include
* **`-f`** - Include

| Base10 Value | Flag          | Description                                                     |
|--------------|---------------|-----------------------------------------------------------------|
| 1            | PAIRED        | paired-end (or multiple-segment) sequencing technology          |
| 2            | PROPER_PAIR   | each segment properly aligned according to the aligner          |
| 4            | UNMAP         | segment unmapped                                                |
| 8            | MUNMAP        | next segment in the template unmapped                           |
| 16           | REVERSE       | SEQ is reverse complemented                                     |
| 32           | MREVERSE      | SEQ of the next segment in the template is reverse complemented |
| 64           | READ1         | the first segment in the template                               |
| 128          | READ2         | the last segment in the template                                |
| 256          | SECONDARY     | secondary alignment                                             |
| 512          | QCFAIL        | not passing quality controls                                    |
| 1024         | DUP           | PCR or optical duplicate                                        |
| 2048         | SUPPLEMENTARY | supplementary alignment                                         |

[REF](http://www.htslib.org/doc/1.11/samtools-flags.html)

#### Paired reads mapped to different scaffolds (`-F 14`)
* Filter out
  * 0x2 each segment properly aligned according to the aligner
  * 0x4 segment unmapped
  * 0x8 next segment in the template unmapped

```
BAM=...
$ samtools view -F 14 ${BAM}
```

### Total count of reads in paried-end BAM - paired vs. unpaired

**Total Reads = Reads_Properly_Paired + Unpaired_Reads**
* Need to consider paired vs. unpaired by filtering on [PROPER_PAIR](http://www.htslib.org/doc/1.11/samtools-flags.html)-flag

```
# Get Total Alignments of Properly-Paired Reads for chromosome X
#    - Note - filter on 3rd column, RNAME, to only grab chrX
$ samtools view -f 0x2 ${BAM} | cut -f3 | grep "chrX" | wc -l
12863955
# Get Total BAM alignments
$ samtools view ${BAM} | cut -f3 | grep "chrX" | wc -l
13246578

# Paired Reads (P) 		= 2 * 12863955		= 25,727,910
# Unpaired Reads (U) 		= 13246578 - 12863955 	= 382,623
# Total Reads (T) = P + U 	= 25,727,910 + 382,623	= 26,110,533
```

## GZIP'd Files can be concatenated
```
$ ls -1
chrX.r1.fq.gz
chrY.r1.fq.gz

$ cat chrX.r1.fq.gz chrY.r1.fq.gz > sex.r1.fq.gz

$ gunzip sex.r1.fq.gz

$ ls *.fq
sex.r1.fq
```
