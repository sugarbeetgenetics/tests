#######UNFINISHED

#This script is to find translocation in TR520 mapping against RefBeetv1.2

import progressbar
import argparse
import os
import numpy as np
import pandas as pd
from datetime import datetime

parser = argparse.ArgumentParser(description='')
parser.add_argument('-mf', '--file', type=str, metavar='', required=True, help='Samtools depth file')
#parser.add_argument('-w', '--window', type=int, metavar='', required=, help='Window size.')
#parser.add_argument('-s', '--step', type=int, metavar='', required=True, help='Window step size.')
args = parser.parse_args()

depth_file = args.file
#win_size = args.window
#win_step = args.step
avg_depth = []
scaffold_ids = []

print("file read!")
print(datetime.now())

with open(depth_file) as d_f:
    read_d_f = pd.read_csv(d_f, sep='\t', header=None, index_col=False)
    data_raw = np.array(read_d_f)
    data = data_raw[data_raw[:, 0].argsort()]
    print(datetime.now())
    print("read csv!")
    #win_start = 0
    #win_end = win_step
    scaffolds = np.unique(data[:, 0])
    print(datetime.now())
    print("entering into for loop!")
    for scaffold in progressbar.progressbar(scaffolds):
        print(datetime.now())
        scaffold_ids.append(scaffold)
        one_scaffold = np.array(np.where(data == scaffold))
        avg_depth.append(float(np.average(one_scaffold[:, 2])))
        data = np.delete(data, np.where(data[:, 0] == scaffold), axis=0) #to make searchable data smaller for next
        # round
            #win_start = win_start + win_step
            #win_end = win_end + (win_size-win_step)
        print("for loop done!")


print(datetime.now())
cwd = os.getcwd()

data = np.column_stack((scaffold_ids, avg_depth))

df = pd.DataFrame(data)

with open(str(cwd+'/'+'svg_depth_per_scaffold.txt'), 'w') as out_file:
    df.to_csv(out_file, sep="\t", header=False, index=False)

print(datetime.now())
