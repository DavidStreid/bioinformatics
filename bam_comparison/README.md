# BAM comparison
Compares the numerical fields of two input BAM files created from different aligners (e.g. compares one aligner's MAPQ score to another on the same input FASTQ files)

Tools:
* [BamUtil](https://genome.sph.umich.edu/wiki/BamUtil:_diff)

## Steps

1. **Sort BAM**
```
samtools sort b1.bam -o b1_sorted.bam
samtools sort b2.bam -o b2_sorted.bam 
```

2. **Run bam diff**
```
bam diff --in1 b1_sorted.bam --in2 b2_sorted.bam --mapQual >> bam_diff.csv
```
* Note - `--mapQual` can be replaced w/ other metric to compare

3. **Parse `bam diff` output**
```
python bam_util_to_csv.py bam_diff.out      # outputs bam_differences.csv
``` 

4. **Graph** 

    BAR
    ```
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
    ```
    
    Scatter (v1 v. v2)
    ```
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import colors
    
    df = pd.read_csv('bam_differences.csv')
    df.plot.scatter(x='v1', y='v2')
    
    plt.title("BAM1 Score vs BAM2 Score")
    plt.xlabel("BAM1 Score")
    plt.ylabel("BAM2 Score")
    ```
    
    Scatter (v1 v. diff)
    ```
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import colors
    
    df = pd.read_csv('bam_differences.csv')
    df.plot.scatter(x='v1', y='v2')
    
    plt.title("BAM1 Score vs (BAM1-BAM2)")
    plt.xlabel("BAM1 Score")
    plt.ylabel("(BAM1-BAM2)")
    ```
