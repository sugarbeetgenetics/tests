
from tqdm import tqdm
import argparse
import os
import numpy as np
import pandas as pd
import re

parser = argparse.ArgumentParser(description='')
parser.add_argument('-l1','--list1', type=str, metavar='', required=True, help='')
parser.add_argument('-s1','--sep1', type=str, metavar='', required=False, help='Set Seperator (Default: Comma)', default="\t")
parser.add_argument('-l2','--list2', type=str, metavar='', required=True, help='')
parser.add_argument('-s2','--sep2', type=str, metavar='', required=False, help='Set Seperator (Default: Comma)', default="\t")
parser.add_argument('-k1','--key1', type=str, metavar='', required=True, help='')
parser.add_argument('-k2','--key2', type=str, metavar='', required=True, help='')
parser.add_argument('-out','--out_file', type=str, metavar='', required=True, help='out file name')
args = parser.parse_args()

version = "merge_tables_v1.1"
cwd = os.getcwd()
out_file_name = cwd + "/" + version + "_" + args.out_file + str(".txt")

list_items1 = args.list1
seperator1 = args.sep1
key1_value = args.key1
key2_value = args.key2

list_items1_df = pd.read_csv(list_items1, sep=seperator1, comment="#", header=[0], low_memory=False)

list_items2 = args.list2
seperator2 = args.sep2

list_items2_df = pd.read_csv(list_items2, delimiter="\t", comment="#", header=[0], low_memory=False)

# merge
merged_table  = list_items1_df.merge(list_items2_df, how="left", left_on =  key1_value, right_on = key2_value)

with open(out_file_name, 'w') as out_file:
    merged_table.to_csv(out_file, sep="\t", index=False)
