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

    grouped = df.groupby('v1').mean()
    grouped['v1-v2'].plot(kind='bar')

    plt.title("BAM Value vs. Average Difference")
    plt.xlabel("BAM1 Score")                    # Or 'BAM2 Score'
    plt.ylabel("Avg(BAM1 Score - BAM2 Score)")  # Avg(BAM2 Score - BAM1 Score) 
    plt.xlim([-2, 62])                          # Helpful for plotting MAPQ, which has range [1-60]

    basename = f.split('.')[0]
    plt.savefig("{}.pdf".format(basename))
