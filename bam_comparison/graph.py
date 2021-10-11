import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.pyplot import figure 
import glob

figure(figsize=(16,8), dpi=80)

files = glob.glob('*_bam_differences.csv')

for f in files:
    df = pd.read_csv(f)

    grouped = df.groupby('control_mapq').mean()
    grouped['ct_diff'].plot(kind='bar')

    plt.title("CTRL-BAM MAPQ vs. Avg( CTRL-BAM MAPQ - EXP-BAM MAPQ )")
    plt.xlabel("CTRL-BAM MAPQ")                    # Or 'BAM2 Score'
    plt.ylabel("Avg( CTRL-BAM MAPQ - EXP-BAM MAPQ )")
    plt.xlim([-2, 62])                          # Helpful for plotting MAPQ, which has range [1-60]

    basename = f.split('.')[0]
    plt.savefig("{}.pdf".format(basename))
