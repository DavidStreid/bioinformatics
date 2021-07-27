# SAM comparison
Compares the numerical fields of two input SAM files created from different aligners (e.g. compares one aligner's MAPQ score to another on the same input FASTQ files)

Tools:
* [BamUtil](https://genome.sph.umich.edu/wiki/BamUtil:_diff)

## Steps

1. Sort BAM
```
samtools sort b1.bam -o b1_sorted.bam
samtools sort b2.bam -o b2_sorted.bam 
```

2. Run bam diff
```
bam diff --in1 b1_sorted.bam --in2 b2_sorted.bam --mapQual >> bam_diff.csv
```
* Note - `--mapQual` can be replaced w/ other metric to compare

3. Parse `bam diff` output
```
python bam_util_to_csv.py bam_diff.out      # outputs bam_differences.csv
``` 

4. Graph
```
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

df = pd.read_csv('bam_differences.csv')
grouped = df.groupby('v1').mean()           # Or 'v2'
grouped['v1-v2'].plot(kind='bar')           # Or 'v2-v1'

plt.title("BAM Value vs. Average Difference")
plt.xlabel("BAM1 Score")                    # Or 'BAM2 Score'
plt.ylabel("Avg(BAM1 Score - BAM2 Score)")  # Avg(BAM2 Score - BAM1 Score) 
```
