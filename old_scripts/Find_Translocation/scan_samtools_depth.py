
import progressbar
import argparse
import os
import numpy as np
import pandas as pd
from datetime import datetime

parser = argparse.ArgumentParser(description='')
parser.add_argument('-mf', '--file', type=str, metavar='', required=True, help='Samtools depth file')

args = parser.parse_args()

depth_file = args.file
result = np.array([["#header", 0, 0]])

print("file read!")
print(datetime.now())


with open(depth_file) as d_f:
    data = np.array(d_f.readlines())
    for line in progressbar.progressbar(data):
        line_split = np.array([str(line).split(sep="\t")])
        #print (str(line_split[0, 2]).split(sep="\n")[0])
        line_split[0,2] = str(line_split[0, 2]).split(sep="\n")[0]
        #check for scaffold, if it is present update it, if not then add new row
        if line_split[:, 0] in result[:, 0]:
            scaf_data = result[np.where(result[:, 0] == line_split[0, 0])]
            scaf_depth = int(scaf_data[0, 1]) + int(line_split[0, 2])
            scaf_freq = int(scaf_data[0, 2]) + 1
            result[np.where(result[:, 0] == line_split[0, 0])[0].tolist(), 0] = line_split[0, 0]
            result[np.where(result[:, 0] == line_split[0, 0])[0].tolist(), 1] = scaf_depth
            result[np.where(result[:, 0] == line_split[0, 0])[0].tolist(), 2] = scaf_freq
        else:
            result = np.vstack((result, line_split))


print(datetime.now())
cwd = os.getcwd()

df = pd.DataFrame(result)

with open(str(cwd+'/'+'svg_depth_per_scaffold.txt'), 'w') as out_file:
    df.to_csv(out_file, sep="\t", header=False, index=False)
print(datetime.now())