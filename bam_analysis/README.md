# BAM ANALYSIS
Collection of scripts to analyze a bam

## TASKS

### `samtools` - piping STDOUT to compressed output

* If outputting to CRAM, include the `-C -T $REF` (e.g. `-C -T hg38.fa`)
* Add `-@ $T` to add more threads (e.g. `@ 16`)
 
Output to BAM

```
./cmd.sh | \
  samtools view -b | \
  tee out.bam | \
  samtools index - out.bam.bai
```

Output to CRAM

```
./cmd.sh | \
  samtools view -C -T $REF - | \
  tee out.cram | \
  samtools index -c - out.cram.crai
```

### Modify readgroups
```
ID=$1
SAM=$2
samtools addreplacerg \
  -r "ID:${ID}" ${SAM} # -O CRAM # If CRAM, must already be in CRAM format
```

## SCRIPTS

###  haplotag_count.sh
Counts the number of reads with haplotag identifiers, such as those added by [whatshap](https://whatshap.readthedocs.io/en/latest/guide.html#phasing-represented-by-hp-tag)
```
./haplotag_count.sh ${BAM}
```
