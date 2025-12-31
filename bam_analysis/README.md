# BAM ANALYSIS
Collection of scripts to analyze a bam

## samtools

### `samtools` - piping STDOUT to compressed output

* If outputting to CRAM, include the `-C -T $REF` (e.g. `-C -T hg38.fa`)
* Add `-@ $T` to add more threads (e.g. `@ 16`)
 
Output to BAM

```
OUT=$1
./cmd.sh | \
  samtools view -b | \
  tee ${OUT} | \
  samtools index - ${OUT}.bai
```

Output to CRAM

```
OUT=$1
REF=$2
./cmd.sh | \
  samtools view -C -T $REF - | \
  tee ${OUT} | \
  samtools index -c - ${OUT}.crai
```

### `samtools tview` - view SAM file like IGV

[tview docs](https://www.htslib.org/doc/samtools-tview.html)

```
# e.g. view SAMPLE alignments in SRY gene region
samtools tview -p chrY:2786855 sample.bam
```
![samtools tview](https://samtools.sourceforge.net/images/seq2-156.png)

### `samtools depth` - retrieve depth of a SAM at a specific region

```
$ samtools depth -r chrM:11800-12000 sample.bam
chrM	11800	201
chrM	11801	200
...
chrM	11999	160
chrM	12000	161
```

### `samtools addreplacerg` - Modify readgroups
* `-w` is to overwrite the readgroup if it already exists
* can also add `PL` for the sequencing platform (e.g. `\tPL:ILLUMINA`)

```
ID=$1
SAM=$2
samtools addreplacerg \
  -w -r "ID:${ID}\tLB:Seq\tSM:${ID}" ${SAM} # -O CRAM # If CRAM, must already be in CRAM format
```

* See `./update_rg.sh` - takes new ID and SAM file and writes to BAM (can write to CRAM, just needs manual edit of `SFX`)

### `samtools collate` / `samtools fastq` - BAM -> FASTQ

* `samtools collate` - outputs SAM grouped by read-name
* `samtools fastq` - convert BAM to FASTQ

```
BAM=...
samtools collate ${BAM} -u -@ 8 -O | \
  samtools fastq -1 R1.fq -2 R2.fq -n -0 ambiguous.fq -s singletons.fq

# ambiguous.fq - R1 & R2 flags are both set or unset  (or /dev/null)
# singleton.fq - Unpaired reads                       (or /dev/null)
```

## bedtools

### BAM -> BED

#### `bamtobed` - retrieves **MAPPED** regions covered by BAM

* **IMPORTANT**: `bamtobed` will exclude hard-clipped, soft-clipped, and inserted regions as well as unmapped reads
  * e.g. `CIGAR=59S48M2D25M19S` -> 75 bp region in the bed

```
bedtools bamtobed -i SAMPLE.bam | bedtools merge -i stdin > SAMPLE.merged.bed
```

## SCRIPTS

###  haplotag_count.sh
Counts the number of reads with haplotag identifiers, such as those added by [whatshap](https://whatshap.readthedocs.io/en/latest/guide.html#phasing-represented-by-hp-tag)
```
./haplotag_count.sh ${BAM}
```
