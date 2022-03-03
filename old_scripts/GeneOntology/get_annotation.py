
from tqdm import tqdm
import argparse
import os
import numpy as np
import pandas as pd
import re

from goatools import obo_parser
p = obo_parser.GODag('/media/pflanz/Hs1-2/Scripts/GeneOntology/go-basic.obo') # get go database


version="v1"
cwd = os.getcwd() + "/" + version

parser = argparse.ArgumentParser(description='Search for go annotation')
parser.add_argument('-l','--list', type=str, metavar='', required=True, help='go ids list')

args = parser.parse_args()

list_items = args.list
list_items_df = pd.read_table(list_items, delimiter="\t", comment="#", names=["ids"], low_memory=False)
list_items_df_np = np.array(list_items_df).squeeze()

names = []
#name = []

counter = 1
for i in tqdm(range(0,len(list_items_df_np[:,1]),1)):
    ids = list_items_df_np[i,1]
    names.append(p[str(ids)].name)

new_s = ','.join(names)

with open(str(cwd + '_' + 'go_Annotation.out'), 'w') as out_file:
    out_file.write(new_s)
