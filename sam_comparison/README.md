## THIS WAS AN EDUCATIONAL EXPERIENCE. In the future, just use BamUtil

[BamUtil](https://genome.sph.umich.edu/wiki/BamUtil:_diff)

# Sam comparison
Compares the numerical fields of two input SAM files created from different aligners (e.g. compares one aligner's MAPQ score to another on the same input FASTQ files)

## Notes
* Identifies matching Reads between two SAM files by their `QNAME`. Paired-End reads will have the same template and therefore the same `QNAME` so to identify the correct read to compare, the `FLAG` value is used. The R1 of each `.sam` is identified by a `FLAG < 128` and the R2 is identified by a `FLAG > 127`. See [Sequence Alignment/Map Format Specification](https://samtools.github.io/hts-specs/SAMv1.pdf) for more information on SAM fields.
