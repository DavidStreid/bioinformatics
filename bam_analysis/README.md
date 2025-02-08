# BAM ANALYSIS
Collection of scripts to analyze a bam

## TASKS

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

### Modify readgroups
* `-w` is to overwrite the readgroup if it already exists
* can also add `PL` for the sequencing platform (e.g. `\tPL:ILLUMINA`)

```
ID=$1
SAM=$2
samtools addreplacerg \
  -w -r "ID:${ID}\tLB:Seq\tSM:${ID}" ${SAM} # -O CRAM # If CRAM, must already be in CRAM format
```

## SCRIPTS

###  haplotag_count.sh
Counts the number of reads with haplotag identifiers, such as those added by [whatshap](https://whatshap.readthedocs.io/en/latest/guide.html#phasing-represented-by-hp-tag)
```
./haplotag_count.sh ${BAM}
```
