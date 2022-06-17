# Simple
Simple things 

## Extracting R1/R2 FASTQs from BAMs
* Rembmer to sort BAM by name so that paired reads aren't missed

```
echo "Sorting SAM=${sorted_sam}"
samtools sort -n -o ${sorted_sam} ${original_sam}
echo "Extracting reads in FQ=${fq1} FQ2=${fq2}"
bedtools bamtofastq -i ${sorted_sam} -fq ${fq1} -fq2 ${fq2}
```

* Verify
```
samtools flagstat ${original_sam}
seqkit stats ${fq1}
seqkit stats ${fq2}
```

