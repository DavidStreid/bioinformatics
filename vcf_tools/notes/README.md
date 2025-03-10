# Notes

## FORMAT

### PL Ordering
Phred Genotype Likelihood scores are ordered by default according to the [VCF Format Specifications](https://samtools.github.io/hts-specs/VCFv4.2.pdf), which provides the following equation

```
F(j/k) = (k*(k+1)/2)+j
where k=allele_2, j=allele_1
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
0       1       2       3       4       5       6       7       8       9       10      11      12      13      14
R/R     R/A     A/A     R/G     A/G     G/G     R/C     A/C     G/C     C/C     R/T     A/T     G/T     C/T     T/T
```

```
R/R     (0 * (0 + 1) / 2) + 0 = 0
R/A     (1 * (1 + 1) / 2) + 0 = 1
A/A     (1 * (1 + 1) / 2) + 1 = 2
R/G     (2 * (2 + 1) / 2) + 0 = 3
A/G     (2 * (2 + 1) / 2) + 1 = 4
G/G     (2 * (2 + 1) / 2) + 2 = 5
R/C     (3 * (3 + 1) / 2) + 0 = 6
A/C     (3 * (3 + 1) / 2) + 1 = 7
G/C     (3 * (3 + 1) / 2) + 2 = 8
C/C     (3 * (3 + 1) / 2) + 3 = 9
R/T     (4 * (4 + 1) / 2) + 0 = 10
A/T     (4 * (4 + 1) / 2) + 1 = 11
G/T     (4 * (4 + 1) / 2) + 2 = 12
C/T     (4 * (4 + 1) / 2) + 3 = 13
T/T     (4 * (4 + 1) / 2) + 4 = 14
```



