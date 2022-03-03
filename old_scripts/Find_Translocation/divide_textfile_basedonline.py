

#This script can divide big csv file for into pieces but none of the piece would divide any scaffold

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


cwd = os.getcwd()

#define funciton to write data to new file in loop


def writer(file_name, data_to_write):
    df = pd.DataFrame(data_to_write)
    with open(str(file_name), 'w') as out_file:
        df.to_csv(out_file, sep="\t", header=True, index=False)


parser = argparse.ArgumentParser(description='Divide big csv file into pieces in a way that none of the piece '
                                             'would divide any scaffold.')
parser.add_argument('-f', '--file', type=str, metavar='', required=True, help='file to be divided')
parser.add_argument('-fs', '--file_size', type=int, metavar='', required=True, help='file size(number of lines)')
args = parser.parse_args()

data_file = args.file
file_size = args.file_size

data_subset_scaf = []
data_subset_base = []
data_subset_depth = []
subset = []
file_num = 0
scaf_change_position = None

with open(data_file) as d_f:
    print("Started processing at %s" %(datetime.now()))
    i = 0
    for line in progressbar.progressbar(d_f):
        data_subset_scaf.append(line.split(sep='\t')[0])
        data_subset_base.append(line.split(sep='\t')[1])
        data_subset_depth.append((line.split(sep='\t')[2]).split(sep='\n')[0])
        subset.append(line.split(sep='\t')[0])
        if i <= file_size:
            if len(subset) >= 2 and subset[-2] == subset[-1]:
                pass
            elif len(subset) >= 2 and subset[-2] != subset[-1]:
                scaf_change_position = i
        elif i > file_size:
            name = cwd+"/subfile"+str(file_num)+".txt"
            #print("\nWriting file! at %s" %(name))
            d_s = pd.DataFrame()
            d_s['Scaffold'] = data_subset_scaf
            d_s['base'] = data_subset_base
            d_s['depth'] = data_subset_depth

            writer(name, d_s[0:scaf_change_position])
            data_subset_scaf = data_subset_scaf[scaf_change_position:]
            data_subset_base = data_subset_base[scaf_change_position:]
            data_subset_depth = data_subset_depth[scaf_change_position:]
            subset = subset[scaf_change_position:]
            file_num += 1
            i = i - scaf_change_position+1
        i +=1
    name = cwd + "/subfile" + str(file_num) + ".txt"
    d_s = pd.DataFrame()
    d_s['Scaffold'] = data_subset_scaf
    d_s['base'] = data_subset_base
    d_s['depth'] = data_subset_depth

    writer(name, d_s[0:scaf_change_position])
