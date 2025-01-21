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
grep -Ev "^#" sample.vcf  | awk '{print $1"\t"$2 - 1"\t"length($4) + $2}' # > sample.bed && bedtools merge -i sample.bed
```

![0/1-based](http://s16.postimg.cc/9ne4syrp1/insertion_or_deletion.jpg)
* Notes - VCF is 1-based and BED is 0-based (See [this post](https://www.biostars.org/p/84686/#84686))

### NOTE on VCF/BED `bedtools intersect`
* `bedtools intersect` will take the length of the REF into account. As in if the BED region position does not intersect `POS`, but does intersect `POS + len(ref)`, then the entry will be output in the intersection. See below -

```
$ tail -1 test.bed test.vcf
==> test.bed <==
chr1	2	3

==> test.vcf <==
chr1	1	v1	TAC	T	.	PASS	.	.	.
$ bedtools intersect -a test.vcf -b test.bed
chr1	1	v1	TAC	T	.	PASS	.	.	.
$ bedtools intersect -a test.vcf -b test.bed
chr1	1	v1	TAC	T	.	PASS	.	.	.
$ bedtools intersect -f 0.33 -a test.vcf -b test.bed
chr1	1	v1	TAC	T	.	PASS	.	.	.
$ bedtools intersect -f 0.34 -a test.vcf -b test.bed
$ 
```


