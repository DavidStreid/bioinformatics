# Bedtools 
Notes for using bedtools with VCFs

## INPUT
```
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	26282344	v_1	ACC	A	.	PASS	.	GT	1/1
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	26282346	v_2	C	G	.	PASS	.	GT	1/1
```

## SUBTRACT
Subtracts from A all features in B that **overlap**, which can be affected by `-f`, `-F`, `-A`, & `-r`

Enforcing overlap percentage w/ `-f`,

* Adding this option **makes subtract more strict**. The overlap needs to meet this proportion in order to remove it from output. Like below where the DEL will cover 100% of the SNV, but not vice-versa.
* `-F` would be the same, except would be as a fraction of B (e.g. A=SNV, B=DEL would output the variant)

```
$ bedtools subtract -f 1 -a del.vcf -b snv.vcf
chr1	26282344	v_1	ACC	A	.	PASS	.	GT	1/1
$ bedtools subtract -f 1 -b del.vcf -a snv.vcf
$ 
```

To make this the most strict, use `-A`, which will not output ANY feature that has any overlap with entries found in B
```
$ bedtools subtract -A -a snv.vcf -b del.vcf
$ bedtools subtract -A -b snv.vcf -a del.vcf
$ 
```

## INTERSECT
You can get the exact complement of `bedtools subtract -A -a $v1 -b $v2` with `bedtools intersect -a $v1 -b $v2` without any options. E.g.
```
$ bedtools intersect -a del.vcf  -b snv.vcf
chr1	26282344	v_1	ACC	A	.	PASS	.	GT	1/1
$ bedtools intersect -b del.vcf  -a snv.vcf
chr1	26282346	v_2	C	G	.	PASS	.	GT	1/1
```

