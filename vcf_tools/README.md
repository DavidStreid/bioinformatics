# VCF Tools
Scripts and notes for analyzing VCF files

## Indexing
* `bgzip`/`tabix` - important for indexing a file for easy viewing
```
bgzip ${VCF}
tabix -p vcf ${VCF}.gz
bcftools view -r chr1:10000-1000000 ${VCF}.gz -o subset_chr1_10000_1000000.vcf.gz
```
