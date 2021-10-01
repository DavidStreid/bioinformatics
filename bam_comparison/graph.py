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
