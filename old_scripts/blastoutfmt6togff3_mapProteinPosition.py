import argparse
import os
#import csv
import re
import numpy as np
import pandas as pd

cwd = os.getcwd()

parser = argparse.ArgumentParser(description='')
parser.add_argument('-b_f','--blast_f', type=str, metavar='', required=True, help='blast file')
parser.add_argument('-r_f','--ref_f', type=str, metavar='', required=True, help='reference file')
args = parser.parse_args()

open_b_f = args.blast_f
open_r_f = args.ref_f

bl_fl = pd.read_table(open_b_f, delimiter="\t", comment="#", names=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])
rf_fl = pd.read_table(open_r_f, delimiter="\t", comment="#", names=["chr","start","end","strand","id"])

bl_fl_np = np.array(bl_fl)
rf_fl_np = np.array(rf_fl)
new_df = np.array(bl_fl_np)

for i in range(0, bl_fl_np.shape[0],1):
    start = int(rf_fl_np[rf_fl_np[:,4]==bl_fl_np[i,0]][:,1])
    end = int(rf_fl_np[rf_fl_np[:,4]==bl_fl_np[i,0]][:,2])
    strand = rf_fl_np[rf_fl_np[:,4]==bl_fl_np[i,0]][:,3]
#    if (strand=="+"):
    new_df[i,0] = "chr6H"
    new_df[i,6] = start + (3*bl_fl_np[i,6])
    new_df[i,7] = end - (3*bl_fl_np[i,7])
    
new_df[:,1]=bl_fl_np[:,1]
new_df[:,2]=bl_fl_np[:,2]
new_df[:,3]=bl_fl_np[:,3]
new_df[:,4]=bl_fl_np[:,4]
new_df[:,5]=bl_fl_np[:,5]
new_df[:,8]=bl_fl_np[:,8]
new_df[:,9]=bl_fl_np[:,9]
new_df[:,10]=bl_fl_np[:,10]
new_df[:,11]=bl_fl_np[:,11]

new_df_pd = pd.DataFrame(new_df, columns = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])

with open(str(cwd + '/' + 'proteins_blastp_put_mapped.txt'), 'w') as out_file:
    new_df_pd.to_csv(out_file, sep="\t", header=True, index=False)