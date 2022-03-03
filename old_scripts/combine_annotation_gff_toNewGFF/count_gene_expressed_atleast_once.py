from tqdm import tqdm
import argparse
import os
#import csv
import re
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


version="v1"
cwd = os.getcwd() + "/" + version

parser = argparse.ArgumentParser(description='id for gene which are expressed atleast once')
parser.add_argument('-fpkm_all','--fpkm_all', type=str, metavar='', required=True, help='fpkm_tracking all')

args = parser.parse_args()

fpkm_all_file = args.fpkm_all
fpkm_all_df = pd.read_csv(fpkm_all_file, delimiter="\t", comment="#", header=0)

total_failed = 0
for i in list(fpkm_all_df['tracking_id'].unique()):
    #print (i)
    all_stages = fpkm_all_df.loc[fpkm_all_df['tracking_id']==i,['tracking_id', 'q1_FPKM','q2_FPKM','q3_FPKM','q4_FPKM','q5_FPKM','q6_FPKM','q7_FPKM','q8_FPKM']]
    all_stages_list = list(all_stages.iloc[0,:])
    #print(all_stages)
    #print(np.array(all_stages_list).reshape(1,9))
    value_count = 0
    for value in all_stages_list[1:]:
        if float(value) >= 1:
            value_count += 1
    if value_count >= 0:
        print(i)
    #print(total_failed)
