import matplotlib.pyplot as plt
import argparse
import numpy as np
import os
import pandas as pd
import math

version="v1"
cwd = os.getcwd() + "/" + version

parser = argparse.ArgumentParser(description='Bar Plot!')
parser.add_argument('-f','--input_file', type=str, metavar='', required=True, help='input file')
parser.add_argument('-o','--out', type=str, metavar='', required=True, help='outfile')
parser.add_argument('-y1max','--y1max', type=float, metavar='', required=True, help='start of break')
parser.add_argument('-y2min','--y2min', type=float, metavar='', required=True, help='end of break')


args = parser.parse_args()
outfile_name = args.out
input_file = args.input_file
break_start = args.y1max
break_end = args.y2min

list_items_df = pd.read_csv(input_file, sep=",", comment="#", low_memory=False)

list_items_df_p1 = list_items_df[list_items_df['genotype']=="SFAR4a1"]
p1_data = list_items_df_p1.loc[:,'expression']
p1_sem = list_items_df_p1.loc[:,'SEM']
list_items_df_p2 = list_items_df[list_items_df['genotype']=="SFAR4a2"]
p2_data = list_items_df_p2.loc[:,'expression']
p2_sem = list_items_df_p2.loc[:,'SEM']

list_items_df_p3 = list_items_df[list_items_df['genotype']=="SFAR4b1"]
p3_data = list_items_df_p3.loc[:,'expression']
p3_sem = list_items_df_p3.loc[:,'SEM']
list_items_df_p4 = list_items_df[list_items_df['genotype']=="SFAR4b2"]
p4_data = list_items_df_p4.loc[:,'expression']
p4_sem = list_items_df_p4.loc[:,'SEM']


# set width of bar
barWidth = 0.25

# Set position of bar on X axis
r1 = np.arange(len(list_items_df_p1))
r2 = [x + barWidth for x in r1]
plt.title(outfile_name)


#### make sublots
f, (ax, ax2) = plt.subplots(2, 1, sharex=True)
# Make the plot
ax.bar(r1, p1_data, yerr = p1_sem, capsize=3, color='#7f6d5f', width=barWidth, edgecolor='white', label='SFAR1a')
ax.bar(r2, p2_data, yerr = p2_sem, capsize=3, color='#557f2d', width=barWidth, edgecolor='white', label='SFAR1c')

# Make the plot
ax2.bar(r1, p1_data, yerr = p1_sem, capsize=3, color='#7f6d5f', width=barWidth, edgecolor='white', label='SFAR1a')
ax2.bar(r2, p2_data, yerr = p2_sem, capsize=3, color='#557f2d', width=barWidth, edgecolor='white', label='SFAR1c')

ax2.xticks([r + barWidth for r in range(len(list_items_df_p1))], list_items_df['DAP'])

"""
##### Broken axis
"""
# zoom-in / limit the view to different portions of the data
ax.set_ylim(.78, 1.)  # outliers only
ax2.set_ylim(0, )  # most of the data

# hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()


d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal



kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

plt.legend(bbox_to_anchor=(-2, 2), loc='upper right')
plt.savefig(cwd+outfile_name+'.png', dpi=300)
