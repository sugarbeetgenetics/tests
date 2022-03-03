

import os
import numpy as np
import pandas as pd
import progressbar
import argparse

parser = argparse.ArgumentParser(description='Create windows for all the chromosomes.')
parser.add_argument('-f', '--file', type=str, metavar='', required=True, help='Chromosome length file')
parser.add_argument('-w', '--window', type=int, metavar='', required=True, help='Window size')
args = parser.parse_args()

chromosome_info_file = args.file
window_size = args.window

def Create_Windows(window_size, chromosome_info):
    start = []
    end = []
    chr = []
    chr.append(chromosome_info[0])
    start.append(0)
    end.append(window_size)
    while (end[-1]+window_size) < int(chromosome_info[1]):
        start.append(end[-1])
        end.append(end[-1]+window_size)
        chr.append(chromosome_info[0])
    if end[-1] < int(chromosome_info[1]):
        start.append(end[-1])
        end.append(chromosome_info[1])
        chr.append(chromosome_info[0])
    return chr, start, end

all_chr = []
all_starts = []
all_ends = []
with open(chromosome_info_file) as c_i_f:
    for chr in progressbar.progressbar(c_i_f):
        chr_split = chr.split(sep="\n")[0]
        chr_split = chr_split.split(sep="\t")
        chr, start, end = Create_Windows(window_size, chr_split)
        all_starts = all_starts + start
        all_ends = all_ends + end
        all_chr = all_chr + chr

#print(all_starts)

cwd = os.getcwd()
col0 = np.array(all_chr)
col1 = np.array(all_starts)
col2 = np.array(all_ends)

data = np.column_stack((col0, col1, col2))
#print(data[0:5])
df=pd.DataFrame(data)

with open(str(cwd+'/'+'windows.txt'), 'w') as out_file:
    df.to_csv(out_file,sep="\t",header=False,index=False)