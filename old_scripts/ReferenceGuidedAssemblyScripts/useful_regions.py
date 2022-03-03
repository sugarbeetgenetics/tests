#This script generate blocks based on base depth data

#process row by row. store start and stop of a block in a list until reach the super block threshold. once it
#reaches the threshold, store block start i.e. first entry of block_start list and block_end i.e. second last entry of
#block_end list. in the same "if block" reset block_start and block_end list but include last two entries in each in
#order to get overlapping superblocks.

import progressbar
import argparse
import os
#import csv
import re
import numpy as np
import pandas as pd
from datetime import datetime

parser = argparse.ArgumentParser(description='Give genome regions with continuous minimum required coverage')
parser.add_argument('-f', '--file', type=str, metavar='', required=True, help='blocks file')
parser.add_argument('-c', '--coverage', type=int, metavar='', required=True, help='minimum coverage')
parser.add_argument('-l', '--length', type=int, metavar='', required=True, help='minimum length')
args = parser.parse_args()

depth_file = args.file
cov = args.coverage
length = args.length
all_blocks = np.empty((0, 5))


with open(depth_file) as d_f:
    read_d_f = pd.read_csv(d_f, delimiter='\t', header=None)
    file_np = np.array(read_d_f)
    #print(file_np.shape)

    temp_row = file_np[np.where(file_np[:, 3]/(file_np[:, 2]-file_np[:, 1]) >= cov), :][0]
    print(temp_row.shape)
    print(temp_row[0:5, :])
    avg_depth = []

    for row in temp_row:
        avg_depth.append(row[3]/(row[2]-row[1]))
    #print(np.array(key_row))
    key_row_np = np.array(avg_depth)
    key_row_np.shape = (len(avg_depth), 1)
    print(key_row_np)
    update_temp_row = np.concatenate((temp_row, key_row_np), axis=1)
    pd_sorted_temp_row = pd.DataFrame(update_temp_row, columns=['scaffold', 'start', 'end', 'size', 'avg_depth'])
    np_sorted_temp_row = pd_sorted_temp_row.sort_values(['scaffold', 'start'], ascending=[True, True])
    np_sorted_temp_row = np.array(np_sorted_temp_row)
    print(np_sorted_temp_row.shape)

    previous_row = np.array([])
    block_s = np.empty((0, np_sorted_temp_row.shape[1]))

    i = 0
    for row in progressbar.progressbar(np_sorted_temp_row):

        if i > 0:

            if row[0] == previous_row[0] and row[1] == previous_row[2]:
                print('row 1 is %s and previous_row 2 is %s' %(row[1], previous_row[2]))
                block_s = np.append(block_s, [previous_row], axis=0)
                print('line 60')
                print(block_s)

            elif row[0] != previous_row[0] or row[1] != previous_row[2]:
                block_s = np.append(block_s,  [previous_row], axis=0)
                print('line 64')
                print(block_s)
                print(block_s.shape)
                if block_s.shape[0] >= 2:
                    #print(block_s)
                    if block_s[-1, 2] - block_s[0, 1] >= length:
                        all_blocks = np.vstack((all_blocks, block_s))

                block_s = np.array(())
                block_s.shape = (0, np_sorted_temp_row.shape[1])
        #print(row)
        previous_row = row
        i +=1
print("Processing finished at %s" % (datetime.now()))

#print(all_blocks)
print("Printing finished at %s" % (datetime.now()))

cwd = os.getcwd()
all_blocks = np.array(all_blocks)
#print ("%s %s %s %s if" %(len(spblock_sf),len(spblock_s),len(spblock_e), len(spblock_length)))
#data = np.array((all_blocks))
#print(data[0:5])
df=pd.DataFrame(all_blocks)

with open(str(cwd+'/'+'useful_regions.txt'), 'w') as out_file:
    df.to_csv(out_file, sep="\t", header=False, index=False)
