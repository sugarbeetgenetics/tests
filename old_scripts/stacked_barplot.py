import matplotlib.pyplot as plt
import argparse
import numpy as np
import os
import pandas as pd
import math

version="v1"
cwd = os.getcwd() + "/" + version

parser = argparse.ArgumentParser(description='Stacked Bar Plot!')
parser.add_argument('-f','--input_file', type=str, metavar='', required=True, help='input data file')
parser.add_argument('-o','--out', type=str, metavar='', required=False, help='outfile')
parser.add_argument('-hd','--header_file', type=str, metavar='', required=True, help='header file')
parser.add_argument('-chr','--chr', type=str, metavar='', required=True, help='chromosome')
parser.add_argument('-start','--start', type=int, metavar='', required=True, help='start')
parser.add_argument('-stop','--stop', type=int, metavar='', required=True, help='stop')


args = parser.parse_args()
outfile_name = args.out
input_file = args.input_file
header_file = args.header_file
chr = args.chr
start = args.start
stop = args.stop

header_df = pd.read_csv(header_file, sep="\t", comment="#", low_memory=False, header=None)
#print(list(header_df.iloc[:,0]))
# sort header for winter and spring
header_df_sort = header_df.sort_values(header_df.columns[1])

list_items_df = pd.read_csv(input_file, sep="\t", comment="#", low_memory=False, names=list(header_df.iloc[:,0]))

total_accesions = len(list(list_items_df.columns))

#print(total_accesions)


# set width of bar
barWidth = 0.15

plt.title(outfile_name)


######### subset for locus #########
idx_all = []
for idx, item in enumerate(zip(list_items_df.iloc[:,0]==chr,list_items_df.iloc[:,1]>=start, list_items_df.iloc[:,1]<=stop)):
    if item[0]==item[1]==item[2]:
        idx_all.append(idx)
#print(idx_all)

list_items_df_subset = list_items_df.iloc[idx_all,:]


####################################

'''
# use dictionary to create multiple plots in loop
x = list_items_df.iloc[0:500,1]
y = list_items_df.iloc[0:500,2:6]

def loop_plot(x,plots):
    figs={}
    axs={}
    for idx,plot in enumerate(list(plots.columns)):
        y=plots.loc[:,plot]
        figs[idx]=plt.figure()
        axs[idx]=figs[idx].add_subplot(111)
        axs[idx].bar(x, y, width=1)
        axs[idx].set_ylabel = plot
        if idx>0:
            axs[idx].spines['top'].set_visible(False)
    return figs, axs

figs, axs = loop_plot(x,y)
'''
x = list(list_items_df_subset.iloc[:,1])
y = list_items_df_subset.iloc[:,2:]

i=0
fig, (axs1,axs2) = plt.subplots(2)
for plot in list(header_df_sort.iloc[2:,0]):

    if header_df_sort.iloc[2+i,1] == 'winter':
        col = "#333333"
        y1=y.loc[:,plot]
        axs1.plot(x, y1, color=col, markevery=10) #, color=col
        axs1.fill_between(x,y1,0, color=col, alpha =0.5) #, color = col
    elif header_df_sort.iloc[2+i,1] == 'spring':
        col = "chocolate"
        y1=y.loc[:,plot]
        axs2.plot(x, y1, color=col, markevery=10) #, color=col
        #axs2.xaxis.set_ticks(np.arange(start, stop, 50))
        #axs2.set_xticklabels(np.arange(start, stop, 50),rotation=45)
        axs2.fill_between(x,y1,0, color=col, alpha =0.5) #, color = col
    i+=1

axs1.spines['bottom'].set_visible(False)
axs1.spines['top'].set_visible(False)
axs1.spines['right'].set_visible(False)
axs1.get_xaxis().set_visible(False)

axs2.spines['top'].set_visible(False)
axs2.spines['top'].set_visible(False)
axs2.spines['right'].set_visible(False)

#axs2.get_xaxis().set_visible(False)
plt.savefig(cwd+outfile_name+'.png', dpi=400)
plt.clf()