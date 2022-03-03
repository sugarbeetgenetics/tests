#based on scan_samtools_depth_v2.py

#######FINISHED

#This script is to give depth per window

import progressbar
import argparse
import os
import numpy as np
import pandas as pd
from datetime import datetime


result = np.array([["#header", 0, 0]])

def split(container, count):
    #function splitting a container into equal sizes
    return [container[_i::count] for _i in range(count)]

#print("Processing Started at %s" % (datetime.now()))


parser = argparse.ArgumentParser(description='Generate depth for all the windows.')
parser.add_argument('-ws', '--win_size', type=int, metavar='', required=True, help='Window size')
parser.add_argument('-df', '--depth_file', type=str, metavar='', required=True, help='Depth file')
args = parser.parse_args()

win_size = args.win_size
depth_file = args.depth_file
file_name = depth_file + "_Windowed_Depth_Summary.txt"
win_length = []
depth = []
scaffold = []
start = []
end = []

with open(depth_file) as d_f:
    #d_f_data = d_f.readlines()
    d_f_data_1 = np.array(pd.read_csv(d_f, delimiter='\t', header=0))
    print("file read at %s" % (datetime.now()))
#read the whole file
    # find all the unique scaffolds,
        #subset array for each scaffold
            #perform calculation on each subset array
    #print(d_f_data_1[0,:])
    all_scaf = np.unique(d_f_data_1[:, 0])

    for scaf in progressbar.progressbar(all_scaf):
        #print("%s < %s" %(total, d_f_data_1.shape[0]))
        subset = d_f_data_1[np.where(d_f_data_1[:,0]==scaf)]
        win_start = 0
        win_end = win_size
        while win_end in range(subset.shape[0]):
            scaffold.append(scaf)
            start.append(win_start)
            end.append(win_end)
            depth.append(np.sum(d_f_data_1[win_start:win_end, 2]))
            win_start = win_end
            win_end = win_end + win_size
            if win_end > subset.shape[0]:
                scaffold.append(scaf)
                start.append(win_start)
                end.append(subset.shape[0]+1)
                depth.append(np.sum(d_f_data_1[win_start:subset.shape[0]+1, 2]))
                win_start = 0
                win_end = win_size
                break

print("Calculation Done! %s" % (datetime.now()))
cwd = os.getcwd()
df = pd.DataFrame()
df['Scaffold'] = scaffold
df['start'] = start
df['end'] = end
df['Total_Depth'] = depth

with open(str(file_name), 'w') as out_file:
    df.to_csv(out_file, sep="\t", header=True, index=False)
