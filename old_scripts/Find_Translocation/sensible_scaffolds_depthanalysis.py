#######UNFINISHED

#This script is to find meaningful scaffolds by comparing refbeet and Proc mapping depth on TR520

import progressbar
import argparse
import os
import numpy as np
import pandas as pd
from datetime import datetime


#print("Processing Started at %s" % (datetime.now()))


parser = argparse.ArgumentParser(description='Generate depth for all the windows.')
parser.add_argument('-df1', '--depth_file1', type=str, metavar='', required=True, help='Windowed Depth file 1')
parser.add_argument('-df2', '--depth_file2', type=str, metavar='', required=True, help='Windowed Depth file 2')
args = parser.parse_args()

depth_file_1 = args.depth_file1
depth_file_2 = args.depth_file2

out_file_name = depth_file_1.split(sep='.txt')[0] + "_Windowed_Depth_Summary_Combined.txt"
win_length = []
depth1 = []
depth2 = []
scaffold = []
start = []
end = []

with open(depth_file_1) as d_f1:
    with open(depth_file_2) as d_f2:
        d_f1_data = np.array(pd.read_csv(d_f1, delimiter='\t', header=None))
        d_f2_data = np.array(pd.read_csv(d_f2, delimiter='\t', header=None))
    #d_f_data = d_f.readlines()
        print("files read at %s" % (datetime.now()))
        print(d_f1_data[:,0])

        all_scaf_df1 = np.unique(d_f1_data[:, 0])
        all_scaf_df2 = np.unique(d_f2_data[:, 0])
        all_scaf = max(all_scaf_df1.shape[0], all_scaf_df2.shape[0])
        shape = d_f1_data.shape[0]
        col_ID_1 = []
        col_ID_2 = []
        i = 0
        for i in range(d_f1_data.shape[0]):
            col_ID_1.append(str(str(d_f1_data[i,0])+"_"+str(d_f1_data[i,1])+"_"+str(d_f1_data[i,2])))
            #print(str(str(d_f1_data[i,0])+"_"+str(d_f1_data[i,1])+"_"+str(d_f1_data[i,2])))
            i +=1
        i = 0
        for i in range(d_f2_data.shape[0]):
            col_ID_2.append(str(d_f2_data[i,0])+"_"+str(d_f2_data[i,1])+"_"+str(d_f2_data[i,2]))
            i +=1

        col_ID_1 = np.array(col_ID_1)
        col_ID_2 = np.array(col_ID_2)

        for line in progressbar.progressbar(col_ID_1):
            #print("%s < %s" %(total, d_f_data_1.shape[0]))
            if np.where(col_ID_2==line):
                scaffold.append(d_f2_data[np.where(col_ID_2 == line)[0], 0][0])
                start.append(d_f2_data[np.where(col_ID_2 == line)[0], 1][0])
                end.append(d_f2_data[np.where(col_ID_2 == line)[0], 2][0])
                depth1.append(d_f2_data[np.where(col_ID_2==line)[0], 3][0])
                depth2.append(d_f1_data[np.where(col_ID_2 == line)[0], 3][0])


print("Calculation Done! %s" % (datetime.now()))

df = pd.DataFrame()
df['Scaffold'] = scaffold
df['start'] = start
df['end'] = end
df['Depth1'] = depth1
df['Depth2'] = depth2

with open(str(out_file_name), 'w') as out_file:
    df.to_csv(out_file, sep="\t", header=True, index=False)
