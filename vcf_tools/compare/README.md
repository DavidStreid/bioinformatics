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
