# Create Liftover
Simple script that takes source and target reference and creates a liftover chain follow.   

## Overview
The source reference is treated like a reference and the target reference is broken down into its scaffolds and each
is aligned via blat to the source.

## Usage
#### Arguments `(-arg)`
* `-s`: string, source reference FASTA to use 
  > Eg: `-s GRCh37.fa`
* `-t`: string, target reference FASTA to use
  > Eg: `-t hs37d5.fa`
* `-l`: flag (optional), run locally (default: lsf)

Examples:
```
# Submits jobs locally
./create_liftover.sh -s grch37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -t hs37d5/hs37d5.fasta -l

# Submist jobs to LSF
./create_liftover.sh -s grch37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -t hs37d5/hs37d5.fasta
```

## Tools Used

`faSplit`: Splits input reference file into separate files by their scaffold names (e.g. `chr1.fa`, `chr2.fa`, ...)

`faToTwoBit`: Converts a FASTQ file into its twoBit representation, which is an efficient format for storing genomic 
sequences optimized for extremely fast read speeds [REF](https://biojulia.net/BioSequences.jl/v1.0/io/twobit.html)

`twoBitInfo`: Gets info about the 2bit file. For this, the sequence lengths of the split FASTA records are needed to
create the alignment nets

`blat`: Alignment tool used to align an input `.fa` sequence to a reference `.2bit` file

`liftUp`: Map `.lft` of target reference to a new re-mapped `.psl` file using the source references `.psl`

### TODO
axtChain
chainMergeSort
chainNet
netChainSubset
```

## Reference
Steps outlined in [UCSC's "Minimal Steps For LiftOver"](http://genomewiki.ucsc.edu/index.php/Minimal_Steps_For_LiftOver)

