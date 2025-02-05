# BED analysis

Sorting a BED file
```
cat $BED | sort -k1,1V -k2,2n
```

* `-k1,1V` - First (`-k1`) sort by the first column (`,1`) and use version sorting (`-V`) to sort the numbers naturally (e.g. `chr2` comes before `chr10`)

* `-k2,2n` - Second (`k2`) sort by the second column (`,2`) and use numeric sorting (`-n`)


## SUBTRACT
Subtracts from A all features in B that **overlap**, which can be affected by `-f`, `-F`, `-A`, & `-r`

Enforcing overlap percentage w/ `-f`,

* Adding this option **makes subtract more strict**. The overlap needs to meet this proportion in order to remove it from output. Like below where the DEL will cover 100% of the SNV, but not vice-versa.
* `-F` would be the same, except would be as a fraction of B (e.g. A=SNV, B=DEL would output the variant)

### INPUT
```
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	26282344	v_1	ACC	A	.	PASS	.	GT	1/1
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	26282346	v_2	C	G	.	PASS	.	GT	1/1
```

### COMMAND
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
```


## MERGE
`bedtools merge` will include neighboring variants as overlapping, i.e. Merge if `abs(start - end) <= 1`

### COMMAND
```
$ cat pre_merge.bed
chr10	120761742	120761742
chr10	120761744	120761747
chr10	120761748	120761748
$ bedtools merge -i pre_merge.bed
chr10	120761742	120761742  # NOT MERGED - separation greater than 1
chr10	120761744	120761749  # MERGED     - separation equals 1
```
