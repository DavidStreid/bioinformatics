# Simple
Simple things for analyzing common bioinformatic file formats, i.e. SAM/BAM, FASTQ, VCF, BED 

* [SAM/BAM/CRAM](#sambam)
  * [tee index](#tee-index)
  * [sub-sample Alignment (random)](#sub-sample-alignment-randomly)
  * [sub-sample Alignment (by region)](#sub-sample-alignment-by-region)
  * [Filtering SAM Flags](#filtering-sam-flags--f-f)
    * [Primary Alignments Only](#primary-alignments-only)
    * [Improperly-Aligned Reads](#paired-reads-improperly-aligned--f-14)
  * [Count Reads](#total-count-of-reads-in-paried-end-bam---paired-vs-unpaired)
  * [SAM-to-BAM](#sam-to-analysis-ready-bam)
  * [Extracting Reads](#extracting-r1r2-fastqs-from-bams)
  * [BAM-to-CRAM](#bam-to-cram)
  * [CRAM-to-BAM](#cram-to-bam)
  * [CRAM/BAM-to-FASTQ](#cram-or-bam-to-fastq)
* [FASTA/FASTQ](#fasta-and-fastq)
  * [Retrieving FASTA region](#retrieving-fasta-region)
  * [Extracting Components](#extracting-specific-components)
* [VCF](#vcf)
  * [Variant Call Region](#variant-call-region)
  * [Find overlapping variants within a window](#find-overlapping-variants)
  * [Sort/Compress/Index](#sort-compress--index-vcf)
* [BED](#bed)
  * [Get Number of Bases](#get-total-number-of-bases-in-bed-file)
* [OTHER](#other)
  * [Helpful find commands](#helpful-find-commands)
  * [Get Seq Index](#fold---find-aanucleotide-at-position)
  * [Concatenate .gz](#gzipd-files-can-be-concatenated)
  * [Estimate Alignment Time](#estimate-alignment-time---bwa-mem-will-output-the-number-of-processed-reads)
  * [Save time Output to File](#time-command-and-write-to-file)
* [References](#references)

## SAM/BAM
### tee index
`tee` is a helpful command that will pipe input to standard out AND file. This is helpful to write a processed BAM to file and index simultaneously.
```
$ sambamba view -s 0.01 -f bam -h --subsampling-seed 123 ${INPUT} | \
  tee ${OUTPUT} | \
  samtools index - ${OUTPUT}.bai
```

### [Sub-sample Alignment randomly](https://www.biostars.org/p/76791/#76791)
```
# Extract 1% of reads from BAM
$ SS=0.01
$ INPUT=sample.bam
$ OUTPUT=$(basename ${INPUT} | sed "s/.bam/_ss${SS}.bam/g")	# e.g. "sample_ss0.01.bam"

# V1: SAMTOOLS
$ samtools view -s 0.01 -b -h ${INPUT} -o ${OUTPUT}
# V2: SAMBABMBA (Same thing, but faster)
$ sambamba view -s 0.01 -f bam -h --subsampling-seed 123 -o ${OUTPUT} ${INPUT}"
```

### Sub-Sample Alignment by region
```
CRAM=sample.cram        # CRAM to subsample
BED=region.bed          # Region to subsample CRAM
REF=hg38.fa		# FASTA (if converting CRAM - not needed for BAMs)

base=$(basename ${CRAM} | cut -d'.' -f1)
ss_base=${base}__subSampled

# -b: to convert to BAM
samtools view -L ${BED} -@ 8 -T ${REF} ${CRAM} | tee ${ss_base}.cram | samtools index - ${ss_base}.cram.crai

# If input is BAM - much faster
# sambamba view -f bam -h -L ${BED} -t 8 ${BAM} | tee ${ss_base}.bam | samtools index - ${ss_base}.bam.bai 
```

### Filtering SAM Flags (`-F`/`-f`)
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

#### Primary Alignments Only
```
# -F: Exclude FLAG, i.e. output no reads with these FLAGs
#   not primary alignment (0x100)
#   supplementary alignment (0x800)
samtools view -b -F 2304 ${BAM}

# sambamba view -F "(not secondary_alignment) and (not supplementary)" -f bam -h -t 8 ${BAM}
```

#### Paired reads improperly aligned (`-F 14`)
* Filter out
  * 0x2 each segment properly aligned according to the aligner
  * 0x4 segment unmapped
  * 0x8 next segment in the template unmapped

```
BAM=...
$ samtools view -F 14 ${BAM}
```

### Total count of reads in paried-end BAM - paired vs. unpaired
* From [SAMTools, Dave Tang](https://davetang.org/wiki/tiki-index.php?page=SAMTools#Fastest_way_to_count_number_of_reads)
```
#number of reads
samtools idxstats in.bam | awk '{s+=$3+$4} END {print s}'
#number of mapped reads
samtools idxstats in.bam | awk '{s+=$3} END {print s}'
```

#### Total Reads = Reads_Properly_Paired + Unpaired_Reads**
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

### SAM to analysis-ready BAM
* Convert to BAM, sort, & index
```
INPUT_SAM=...
BAM=...         # $(echo ${INPUT_SAM} | sed 's/.sam/.bam/g')    # Get BAM output name
SORT_BAM=...    # $(echo ${BAM} | sed 's/.bam/_sorted.bam/g')   # Get sorted-BAM output name
samtools view -h -b ${INPUT_SAM} > ${BAM}
samtools sort -o ${SORT_BAM} -O BAM
samtools index ${SORT_BAM}
```

### Extracting R1/R2 FASTQs from BAMs
* Rembmer to sort BAM by name so that paired reads aren't missed
* Paired-Read Notes
  * Secondary-alignment mates may be missing in the BAM. To avoid the error below, remove secondary/supplemental alignments, `-F 2304`, i.e. Filter out (`-F`) secondary & supplemntary (`256` + `2048`) reads prior to extracting. 
    ```
    *****WARNING: Query <RGID> is marked as paired, but it's mate does not occur next to it in your BAM file.  Skipping
    ```
  * Will usually end up with the following
    * Paired: Mates present in BAM
    * [Properly-Paired](https://www.biostars.org/p/8318/): 
    
      1. Both mates mapped
      2. Mates map within a reasonable distance of each other (e.g. Different chromosome mapping isn't properly-paired)
      3. Orientation of R1 & R2 is correct, i.e. Forward strand is 5'->3' upstream of reverse strand that is 3'->5'

**1) Extract name-sorted, primary alignments**
```
# OPTION 1) COLLATE
samtools collate ${original_sam} > ${sorted_sam}

# OPTION 2) SORT
echo "Sorting SAM=${sorted_sam}"
samtools view -F 2304 ${original_sam} > primary_alignment.sam
samtools sort -n -o ${sorted_sam} primary_alignment.sam
```

**2) Extract FQ files**
```
echo "Extracting reads in FQ=${fq1} FQ2=${fq2}"
bedtools bamtofastq -i ${sorted_sam} -fq ${fq1} -fq2 ${fq2}
```

**3) Verify**
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

### BAM-to-CRAM 
```
# BAM -> CRAM
$ samtools view -@ 8 -T ${FA_REF} -C -o ${cram} ${bam}

# BAM (sorted) -> CRAM/CRAI
$ samtools view -@ 8 -T ${FA_REF} -C ${bam} | tee ${cram} | samtools index -c - ${cram}.crai
```
* [REF](https://davetang.org/muse/2014/09/26/bam-to-cram/)
* NOTE - Don't specify threads for `samtools index`. It is less compute-intensive than view so best to not take resources from the `view` command. `samtools` 1.11+ defaults to number of cores available and earlier defaults to `-@ 1`
* **[WARNING] Don't use `sambamba`** - This will work, but while `sambamba` usually beats `samtools` in speed on shared functionality, `sambamba` speedup is minimal or nonexistant for CRAMs.
  * This is because `sambamba` wraps the C-language `htslib` library used by `samtools` and doesn't utilize the D-language `sambamba` was built on, which gives it its parallel-processing speed

### CRAM-to-BAM
```
# CRAM -> BAM
$ samtools view -@ 8 -T ${FA_REF} -b ${cram} -o ${bam}

# CRAM (sorted) -> BAM/BAI
$ samtools view -@ 8 -T ${FA_REF} -b ${cram} | tee ${bam} | samtools index - ${bam}.bai
```
* NOTE - Use fewer threads `-@ #` for `samtools index`, or better yet, don't specify threads at all
* **[WARNING] Don't use `sambamba`** - `sambamba` support for CRAMs was removed after v0.8 due to [processing issues](https://github.com/biod/sambamba/issues/280). If trying to convert CRAM-to-BAM w/ `sambamba` it might fail

### CRAM or BAM-to-FASTQ
#### Single-End
* Writes primary alignment reads to FASTQ (`-F 0X900`)
```
# -o: All reads write to file
$ samtools fastq -F 0x900 -o ${FQ} ${CRAM} --reference ${REF}
```
#### Paired-End
```
# -1: reads marked R1 write to file
# -2: reads marked R2 write to file
$ samtools fastq -F 0x900 -1 r1.fq -2 r2.fq ${BAM}
```

## FASTA and FASTQ
### Retrieving fasta region
[bedtools getfasta](https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html)
```
BED=region.bed
REF=hg38.fa

$ bedtools getfasta -fi ${REF} -bed ${BED}
```

### Extracting Specific Components
* Assuming standard [FQ format](https://en.wikipedia.org/wiki/FASTQ_format#Format) with each read & meta-information contained in 4-line components
```
IDX=1
#   0: QUALITY
#   1: SEQ_ID
#   2: SEQUENCE

awk 'NR % 4 == ${IDX}' ${FQ}
```

## VCF
### Variant Call Region
If the region is known, use `bcftools`- 
* `mpileup` to generate a [pileup summary](https://davetang.org/muse/2015/08/26/samtools-mpileup/)
* `call` to call the variant

e.g. Looking for [RFC1 STR](https://github.com/Illumina/ExpansionHunter/blob/master/variant_catalog/grch38/variant_catalog.json#L527)
```
ref=hg38.fa
chr=chr4
start=39348424
end=39348479
aln_file=sample.bam
bcftools mpileup -f {ref} -r {chr}:{start}:{end} {aln_file} | bcftools call -mv -Ov -o sample.vcf 
```

### Find overlapping variants
```
f1=${FN}
f2=${IN}
OVERLAP=5
echo "OVERLAP=${OVERLAP}"
grep -v "^#" ${f1} | awk '{print $1"\t"$2"\t"$2}' > f1.bed
grep -v "^#" ${f2} | awk '{print $1"\t"$2"\t"$2}' > f2.bed

OUT="overlap.bed"
bedtools window -a f1.bed -b f2.bed -w ${OVERLAP} > ${OUT}

bedtools intersect -a ${f1} -b ${OUT} > $(basename ${f1} | sed "s/.vcf/_overlap${OVERLAP}.vcf/g")
bedtools intersect -a ${f2} -b ${OUT} > $(basename ${f2} | sed "s/.vcf/_overlap${OVERLAP}.vcf/g")
```

### Sort, compress, & index VCF
1. Sort (Note - sorts according to the Sequence Dictionary in the VCF's headers)
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

## BED
### Get Total Number of bases in BED file
```
cat file.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
```
* REF - [BioStars](https://www.biostars.org/p/68283/#68292)


## OTHER
### Helpful find commands
**Find all `.bam` & `.bai` files**
```
$ find ./ -type f \( -iname \*.bam -o -iname \*.bai \)
```
**`tar` find output to file**
```
$ find ./ -type f \( -iname \*.vcf -o -iname \*.fastq \) | tar -zcvf secondary_analysis.tar.gz  -T -
```

### Fold - Find AA/Nucleotide at position
e.g. What amino acid is at position 16 of a transcript
```
TR="MITFLPIIFSSLVVVTFVIGNFANGFIALVNSIE"
IDX=16
$ echo "${TR}" | fold -w 1 | cat -n | grep "\s${IDX}\s"
    16	T
```

### GZIP'd Files can be concatenated
```
$ ls -1
chrX.r1.fq.gz
chrY.r1.fq.gz

$ cat chrX.r1.fq.gz chrY.r1.fq.gz > sex.r1.fq.gz

$ gunzip sex.r1.fq.gz

$ ls *.fq
sex.r1.fq
```

### Estimate alignment time - BWA MEM will output the number of processed reads 

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

### Time command AND write to file
Wrap the `time <CMD>` in `()`
```
$ (time ls -l) > output.txt 2>&1
$ echo output.txt
total 20
-rw-r--r--. 1 me grp 284 Oct 11 18:59 output.txt

real    0m0.005s
user    0m0.001s
sys     0m0.004s
```


## References
* [Common Samtools Usage Examples](https://davetang.org/wiki/tiki-index.php?page=SAMTools#Extracting_SAM_entries_mapping_to_a_specific_loci) 
