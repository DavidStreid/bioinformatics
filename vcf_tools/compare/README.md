# Compare VCF files
## By position
```
./compare.sh sample_before.vcf sample_refactor.vcf
vcf1=sample_before.vcf
vcf2=sample_refactor.vcf
vcf1_same=61
vcf2_same=61
vcf1_only=0
vcf2_only=0
jaccard=1       # Samples are exactly the same on position
```

## By position and annotation
* Right now only for VCFs that represent the same sample and have that represented in their headers, specifically the SAMPLE columns
```
$ python compare_vcf_by_annotations.py sample_before.vcf sample_new_annotations.vcf
SAMPLES=['sample1']
SAMPLES=['sample1']
VCF file diffs=84
```

## Note
A valid VCF will specify the type & columns, i.e. add something like this
```
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
```

```
$ ./get_merged_vcf_headers.sh s1.vcf s2.vcf
##fileformat=VCFv4.1
##source=VCF_MERGE
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
```
