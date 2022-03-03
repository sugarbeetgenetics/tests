
import numpy as np
import pandas as pd
import math
import argparse
import os
import sys
from tqdm import tqdm

version="V1"
cwd = os.getcwd()

parser = argparse.ArgumentParser(description='Generate average of a column based on labels from another column.')
parser.add_argument('-i', '--input_data', type=str, metavar='', required=True, help='input data file')
parser.add_argument('-o', '--out_file', type=str, metavar='', required=True, help='out file name')

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

DataFile =  args.input_data

out_file_name = cwd + "/" + version + args.out_file + str(".summary")

line_items_df = pd.read_csv(str(DataFile), sep="\t", comment="#", low_memory=False, names=["name", "bin", "d1", "d2"])
line_items_df_np = np.array(line_items_df)

_, idx = np.unique(line_items_df_np[:,0], return_index=True)
col0_uniq = line_items_df_np[np.sort(idx),0]

all_id_name = []
all_bin_total = []
all_d1_total = []
all_average = []

for col0 in tqdm(col0_uniq):
    subset_data = line_items_df_np[line_items_df_np[:,0]==col0,:]
    bin_total = len(subset_data[:,0])
    d1_total = np.sum(subset_data[:,3])
    
    average = d1_total / len(subset_data[:,0])

    all_id_name.append(col0)
    all_bin_total.append(bin_total)
    all_d1_total.append(d1_total)
    all_average.append(average)

combined_data = pd.DataFrame(list(zip(all_id_name, all_bin_total, all_d1_total, all_average)), columns=['ids', 'total_bin', 'd1_total', 'average'])

with open(str(out_file_name), 'w') as out_file:
    combined_data.to_csv(out_file, sep="\t", index=False)
