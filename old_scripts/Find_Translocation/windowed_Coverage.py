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
parser.add_argument('-wf', '--win_file', type=str, metavar='', required=True, help='Window file')
parser.add_argument('-df', '--depth_file', type=str, metavar='', required=True, help='Depth file')
args = parser.parse_args()

window_file = args.win_file
depth_file = args.depth_file
file_name = window_file + "_Windowed_Depth_Summary.txt"

with open(window_file) as c_i_f:
    data = c_i_f.readlines()
    print("file read at %s" % (datetime.now()))


with open(depth_file) as d_f:
    d_f_data = np.array([d_f.readlines()])
    d_f_data_1 = np.array([l.split(',') for l in d_f_data])
    d_f_data_2 = np.array([l.split('\t') for l in d_f_data_1])
    d_f_data_3 = np.array([l[3].split('\n')[0] for l in d_f_data_2])

    for line in progressbar.progressbar(data):
        line_split = np.empty((1, 3))
        if line == "\n":
            pass
        else:
            line_split = np.array([str(line).split(sep="\t")])
            line_split[0, 2] = str(line_split[0, 2]).split(sep="\n")[0]
            d_f_data_scaf = d_f_data_3[np.where(d_f_data_3[:, 0] == line_split[0, 0])]
            np_start = d_f_data_scaf[(np.where(d_f_data_scaf[:, 1] == line_split[0, 1])[0])[0]]
            np_end = d_f_data_scaf[(np.where(d_f_data_scaf[:, 2] == line_split[0, 2])[0])[0]]
            d_f_subset = d_f_data_scaf[:,np_start:np_end+1]
            depth_per_window = np.sum(int(d_f_subset[:, 2]))
            result = np.array([d_f_data_scaf[0, 0], ])

results = np.array(result)

'''
final_result = np.array([["#header", 0, 0]])

if rank == 0:
    print("final calculation at %s" % (datetime.now()))

    reduced_dim_results = results.reshape((np.shape(results)[0]*np.shape(results)[1], np.shape(results)[2]))
    for row in reduced_dim_results:
        if row[0] in final_result[:, 0]:
            final_result[np.where(final_result[:, 0] == row[0])[0].tolist(), 0] = row[0]
            final_result[np.where(final_result[:, 0] == row[0])[0].tolist(), 1] = int(final_result[np.where(
                final_result[:, 0] == row[0])[0].tolist(), 1]) + int(row[1])
            final_result[np.where(final_result[:, 0] == row[0])[0].tolist(), 2] = int(final_result[np.where(
                final_result[:, 0] == row[0])[0].tolist(), 2]) + int(row[2])
        else:
            final_result = np.vstack((final_result, row))
#    print((final_result))
'''
cwd = os.getcwd()
df = pd.DataFrame(results)
with open(str(file_name), 'w') as out_file:
    df.to_csv(out_file, sep="\t", header=False, index=False)

#    print(datetime.now())
