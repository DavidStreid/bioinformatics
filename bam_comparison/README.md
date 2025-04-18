# BAM comparison
Compares the numerical fields of two input BAM files created from different aligners (e.g. compares one aligner's MAPQ score to another on the same input FASTQ files)

## `bam_comparison_exact.sh`

**USE THIS ONE** - python script will run forever and eventually run out of memory

## `bam_comparison_exact.py`
Compares two BAM files that should be exactly the same - checks all fields and tags. Tags can be out-of-order
* MAJOR RAM restrictions - will load every difference-per-CHROM before releasing memory, e.g. 60GB human BAM -> 200GB RAM usage (for the earlier chroms)
* Outputs a `<CHROM>.columns_only.tsv` & `<CHROM>.values.tsv`, which will have either just the fields, or the fields plus different values, respectively

### RUN
```
$ python bam_comparison_exact.py sample.1.bam sample.2.bam
Processing chrM
	chrM	0
	chrM	250000
	chrM	500000
	chrM	750000
	chrM	1000000
	chrM	1250000
READ_COMPARISON	CHROM=chrM	SHARED=637935	R1_ONLY=0	R2_ONLY=0
Processing chr1
	chr1	1500000
    ...
$ cut -f2,3 bam_comparison.chrM.columns_only.tsv  | grep -v "-" | sort | uniq -c
     12 1	flag
 597917 1	query_qualities
    249 1	tags__MD
    243 1	tags__NM
    206 2	flag,query_qualities
   8981 2	query_qualities
   8981 2	tags__MD
   8903 2	tags__NM
```


## `bam_util_to_csv.py`

Tools:
* [BamUtil](https://genome.sph.umich.edu/wiki/BamUtil:_diff)

### Steps
0. **Install python environment**
```
$ conda create --name bam_compare --file requirements.txt
$ conda activate bam_compare
(bam_compare)$ 
```

1. Run compare
```
CONTROL_BAM=...
TARGET_BAM=...

(bam_compare)$ python bam_util_to_csv.py ${CONTROL_BAM} ${TARGET_BAM}
``` 

2. **Parse `bam diff` output**
```
(bam_compare)$ python bam_util_to_csv.py bam_diff.out      # outputs bam_differences.csv
``` 

3. **Graph** 
```
(bam_compare)$ python graph.py
```

 ___Other Graphs___

BAR

    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import colors
    from matplotlib.pyplot import figure 
    
    figure(figsize=(16,8), dpi=80)
    
    df = pd.read_csv('bam_differences.csv')
    grouped = df.groupby('v1').mean()
    grouped['v1-v2'].plot(kind='bar')
    
    plt.title("BAM Value vs. Average Difference")
    plt.xlabel("BAM1 Score")                    # Or 'BAM2 Score'
    plt.ylabel("Avg(BAM1 Score - BAM2 Score)")  # Avg(BAM2 Score - BAM1 Score) 
    plt.xlim([-2, 62])                          # Helpful for plotting MAPQ, which has range [1-60]
    
    plt.savefig("bam1_v_diff.pdf")
    
SCATTER (v1 v. v2)

    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import colors
    
    df = pd.read_csv('bam_differences.csv')
    df.plot.scatter(x='v1', y='v2')
    
    plt.title("BAM1 Score vs BAM2 Score")
    plt.xlabel("BAM1 Score")
    plt.ylabel("BAM2 Score")
    
SCATTER (v1 v. diff)

    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import colors
    
    df = pd.read_csv('bam_differences.csv')
    df.plot.scatter(x='v1', y='v2')
    
    plt.title("BAM1 Score vs (BAM1-BAM2)")
    plt.xlabel("BAM1 Score")
    plt.ylabel("(BAM1-BAM2)")
