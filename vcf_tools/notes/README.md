# Notes

## FORMAT

### PL Ordering
Phred Genotype Likelihood scores are ordered by default according to the [VCF Format Specifications](https://samtools.github.io/hts-specs/VCFv4.2.pdf), which provides the following equation

```
F(j/k) = (k*(k+1)/2)+j
where k=allele_1, j=allele_2
  e.g.
    biallelic:  AA,AB,BB
    triallelic: AA,AB,BB,AC,BC,CC
```

For instance, given the following gVCF entry with `REF=N` and `ALT=A,G,C,T`,

```
chr1	1000	1	N	A,G,C,T,<*>	GT:PL	2/4:0,1,2,3,4,5,6,7,8,9,10,11,12,13,14
```

The following is how the PL values are assigned.
```
0	1	2	3	4	5	6	7	8	9	10	11	12	13	14
A/A	A/G	G/G	A/C	G/C	C/C	A/T	G/T	C/T	T/T	A/<*>	G/<*>	C/<*>	T/<*>	<*>/<*>
```

```
A/A	(0 *(0 + 1 ) / 2) + 0 = 0
A/G	(1 *(1 + 1 ) / 2) + 0 = 1
A/C	(2 *(2 + 1 ) / 2) + 0 = 3
A/T	(3 *(3 + 1 ) / 2) + 0 = 6
A/<*>	(4 *(4 + 1 ) / 2) + 0 = 10
G/G	(1 *(1 + 1 ) / 2) + 1 = 2
G/C	(2 *(2 + 1 ) / 2) + 1 = 4
G/T	(3 *(3 + 1 ) / 2) + 1 = 7
G/<*>	(4 *(4 + 1 ) / 2) + 1 = 11
C/C	(2 *(2 + 1 ) / 2) + 2 = 5
C/T	(3 *(3 + 1 ) / 2) + 2 = 8
C/<*>	(4 *(4 + 1 ) / 2) + 2 = 12
T/T	(3 *(3 + 1 ) / 2) + 3 = 9
T/<*>	(4 *(4 + 1 ) / 2) + 3 = 13
<*>/<*>	(4 *(4 + 1 ) / 2) + 4 = 14
```



