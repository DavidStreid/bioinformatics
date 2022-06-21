# Simple
Simple things 

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
