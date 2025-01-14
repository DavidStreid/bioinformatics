# VCF Tools
Scripts and notes for analyzing VCF files

## Indexing
* `bgzip`/`tabix` - important for indexing a file for easy viewing
```
bgzip ${VCF}
tabix -p vcf ${VCF}.gz
bcftools view -r chr1:10000-1000000 ${VCF}.gz -o subset_chr1_10000_1000000.vcf.gz
```

## Helpful Scripts
### VCF -> BED

```
grep -Ev "^#" sample.vcf  | awk '{print $1, $2 - 1, length($4) + $2}'
```

![Alt text](http://s16.postimg.cc/9ne4syrp1/insertion_or_deletion.jpg)
* Notes - VCF is 1-based and BED is 0-based (Image from [this post](https://www.biostars.org/p/84686/#84686))


