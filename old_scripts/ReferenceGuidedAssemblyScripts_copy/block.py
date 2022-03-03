#This script generate blocks based on base depth data

import progressbar
import argparse
import os
import csv
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='Create blocks from a mapping file. Block is defined as region with constant coverage.')
parser.add_argument('-f', '--file', type=str, metavar='', required=True, help='samtools depth file')
parser.add_argument('-c', '--coverage', type=int, metavar='', required=True, help='minimum coverage threshold for a base.')
args = parser.parse_args()

depth_file = args.file
block_s = []
block_e = []
block_sf = []
conditional_track = []

with open(depth_file) as d_f:
    read_d_f = csv.reader(d_f, delimiter='\t')
    #i = 0
    in_block = 0
    block_start = []
    block_end = []
    block_scaf = []
    for row in progressbar.progressbar(read_d_f):
        li = str(row).strip()
        if not li.startswith("#"):
            #print (row)
            depth = int(row[2])
            if depth > int(args.coverage):
                if in_block==0:
                    conditional_track.append("a")
                    block_start.append(int(row[1]))
                    block_end.append(int(row[1]))
                    block_scaf.append(row[0])
                    in_block = 1
                    #print("a")
                    #print("%s\t%s" % (block_start, block_end))
                    pass

                elif in_block == 1:
                    conditional_track.append("b")
                    block_end[-1] = block_end[-1] + 1
                    in_block = 1
                    #print("b")
                    #print("%s\t%s" %(block_start,block_end))
                    pass

            elif depth <= int(args.coverage):
                conditional_track.append("c")
                in_block = 0
                #print("c")

            elif block_scaf[-1] != row[0] and conditional_track[-1] != ["c"]:
                conditional_track.append("d")
                in_block = 0


            if conditional_track[-3:] == ["b","b","c"] or conditional_track[-3:] == ["b","b","d"]:
                # to record start and end, every time new block starts
                block_e.append(block_end[-1])
                block_sf.append(block_scaf[-1])
                block_s.append(block_start[-1])
            #print(conditional_track[-3:])

cwd = os.getcwd()

data = np.column_stack((block_sf, block_s, block_e))

df = pd.DataFrame(data)

with open(str(cwd+'/'+'blocks.txt'), 'w') as out_file:
    df.to_csv(out_file,sep="\t",header=False,index=False)