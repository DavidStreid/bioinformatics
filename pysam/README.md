# PYSAM
## QuickStart
```
import pysam

in_aln = '...' # BAM/CRAM/SAM
aln = pysam.AlignmentFile(in_aln, "rb")

# BRCA1 Gene
chr = 'chr17'
start = 4304429543
end=43170245 
reads_in_region = aln.fetch(chr, start, end)
```

## Useful Functions
### [AlignmentFile](https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile)
#### [.fetch()] (https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile.fetch) - Retreive reads in alignment
* Returns iterator over reads
```
reads_in_region = aln.fetch(chr, start, end) 
```

### [AlignedSegment](https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment)
#### [get_blocks](https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_blocks) - Retrieve blocks of alignment
```
blocks = read.get_blocks() # Genomic Coordinates of Alignment
```

#### [get_aligned_pairs](https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_aligned_pairs) - Retrieve per-position alignment of read and reference
* `pysam.aligned_pairs` Format: `[(read_pos, ref_pos, <ref_base>), ...], e.g. [(0, REF_START, <A/G/C/T>), (1, REF_START + 1, <A/G/C/T>), ...]`
* Notes
	* If `with_seq=True` is specified, each tuple will have a third element, `ref_base` with the reference nucleotide
	* Insertions are formatted with `ref_pos` & `ref_base` being `None`
	* Deletions are formatted with `read_pos` being `None`

```
for aln_seg in reads_in_region:
  read.get_aligned_pairs(with_seq=True) # [(int, int, str), ...] <- Has reference nucleotide
  # read.get_aligned_pairs() 		# [(int, int), ...]	 <- No reference nucelotide
```
