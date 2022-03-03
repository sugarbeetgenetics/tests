#This script generate blocks based on base depth data

#process row by row. store start and stop of a block in a list until reach the super block threshold. once it
#reaches the threshold, store block start i.e. first entry of block_start list and block_end i.e. second last entry of
#block_end list. in the same "if block" reset block_start and block_end list but include last two entries in each in
#order to get overlapping superblocks.


import progressbar
import argparse
import os
import csv
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='Create super blocks from block data file. Combine adjacent blocks to '
                                             'generate super blocks upto specified length.')
parser.add_argument('-f', '--file', type=str, metavar='', required=True, help='blocks file')
parser.add_argument('-l', '--length', type=int, metavar='', required=True, help='Maximum length of a super block.')
args = parser.parse_args()

block_file = args.file
spblock_s = []
spblock_e = []
spblock_sf = []
spblock_length = []

with open(block_file) as b_f:
    read_b_f = csv.reader(b_f, delimiter='\t')
    block_start = []
    block_end = []
    block_scaf = []

    for row in progressbar.progressbar(read_b_f):
        block_scaf.append(row[0])
        block_start.append(int(row[1]))
        block_end.append(int(row[2]))



        if block_scaf[0] == block_scaf[-1]:
            sp_block_length = block_end[-1] - block_start[0]
            if sp_block_length > int(args.length):
                spblock_s.append(block_start[0])
                spblock_e.append(block_end[-1])
                spblock_sf.append(block_scaf[0])
                spblock_length.append(sp_block_length)
                #print ("%s %s %s if" %(block_scaf[0],block_start[0],block_end[-1]))
                block_start = []
                block_end = []
                block_scaf = []
                pass
        elif block_scaf[0] != block_scaf[-1]:
            sp_block_length = block_end[-2] - block_start[0]
            spblock_s.append(block_start[0])
            spblock_e.append(block_end[-2])
            spblock_sf.append(block_scaf[0])
            spblock_length.append(sp_block_length)
            #print("%s elif" %(block_start[0]))
            block_start = []
            block_end = []
            block_scaf = []

cwd = os.getcwd()
#print ("%s %s %s %s if" %(len(spblock_sf),len(spblock_s),len(spblock_e), len(spblock_length)))
data = np.column_stack((spblock_sf, spblock_s, spblock_e, spblock_length))
#print(data[0:5])
df=pd.DataFrame(data)

with open(str(cwd+'/'+'superblock.txt'), 'w') as out_file:
    df.to_csv(out_file,sep="\t",header=False,index=False)