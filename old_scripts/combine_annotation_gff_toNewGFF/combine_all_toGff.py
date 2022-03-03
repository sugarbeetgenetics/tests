#Author Avneesh kumar
#Date: 04/01/2019
#This script takes gff and other tables like pfam output or tmhmm output or blast outfmt6 and combine all that info into original gff


import progressbar
import argparse
import os
#import csv
import re
import numpy as np
import pandas as pd
from datetime import datetime

cwd = os.getcwd()

parser = argparse.ArgumentParser(description='Give gff')
parser.add_argument('-gff','--gff', type=str, metavar='', required=True, help='gff')
parser.add_argument('-bp1', '--blastp1', type=str, metavar='', required=True, help='blastp first')
parser.add_argument('-bp2', '--blastp2', type=str, metavar='', required=True, help='blastp second')
parser.add_argument('-pf', '--pfam', type=str, metavar='', required=True, help='pfam out table')
parser.add_argument('-tm', '--tmhmm', type=str, metavar='', required=True, help='tmhmm out table')
parser.add_argument('-sp', '--signalp', type=str, metavar='', required=True, help='signalp out table')
args = parser.parse_args()

gff_file = args.gff
blastp_file_1 = args.blastp1
blastp_file_2 = args.blastp2
pfam_file = args.pfam
tmhmm_file = args.tmhmm
signalp_file = args.signalp

gff_df = pd.read_table(gff_file, delimiter="\t", comment="#", names=["seq_name","source","feature","start","end","score","strand","frame","attribute", "target_gtf"])
blastp_df_1 = pd.read_csv(blastp_file_1, sep="\t", comment="#", names=["qseqid", "sseqid", "pident", "length","e-value"])
blastp_df_2 = pd.read_csv(blastp_file_2, sep="\t", comment="#", names=["qseqid2", "sseqid2", "pident2", "length2","e-value2"])
pfam_df = pd.read_table(pfam_file, delimiter="\t", comment="#", names=["target", "accession1", "query_name"])
tmhmm_df = pd.read_csv(tmhmm_file, sep="\t", comment="#", names=["target_tm", "Topology"])
signalp_df = pd.read_csv(signalp_file, sep="\t", comment="#", names=["target_sp", "SignalP"])

# modify blastp out: get one record for one transcript (best scored).
#convert pd df to np array bcoz it is easier to work with (for me ;) )
blastp_np_1 = np.array(blastp_df_1)
#print(blastp_df.head())
blastp_np_1_filtered=np.empty((0, blastp_np_1.shape[1]))
#print(blastp_np.shape[0])

blastp_np_1_filtered = np.append(blastp_np_1_filtered, [blastp_np_1[0,]], axis=0)

for i in progressbar.progressbar(range(1, blastp_np_1.shape[0],1)):
    if(blastp_np_1[i,0] == blastp_np_1_filtered[-1,0]):
        if(blastp_np_1[i,1] != blastp_np_1_filtered[-1,1]):
            blastp_np_1_filtered[-1,0] = str(blastp_np_1_filtered[-1,0]) + str(',') + str(blastp_np_1[i,0])
    else:
        blastp_np_1_filtered = np.append(blastp_np_1_filtered, [blastp_np_1[i,]], axis=0)

#convert np array to pd df
blastp_df_1_filtered = pd.DataFrame(blastp_np_1_filtered, columns=["qseqid", "sseqid", "pident", "length","e-value"])


with open(str(cwd + '/' + 'blastp_df_1_filtered.out'), 'w') as out_file:
    blastp_df_1_filtered.to_csv(out_file, sep="\t", header=False, index=False)

#convert pd df to np array bcoz it is easier to work with (for me ;) )
blastp_np_2 = np.array(blastp_df_2)
#print(blastp_df.head())
blastp_np_2_filtered=np.empty((0, blastp_np_2.shape[1]))
#print(blastp_np.shape[0])

blastp_np_2_filtered = np.append(blastp_np_2_filtered, [blastp_np_2[0,]], axis=0)

for i in progressbar.progressbar(range(1, blastp_np_2.shape[0],1)):
    if(blastp_np_2[i,0] == blastp_np_2_filtered[-1,0]):
        if(blastp_np_2[i,1] != blastp_np_2_filtered[-1,1]):
            blastp_np_2_filtered[-1,0] = str(blastp_np_2_filtered[-1,0]) + str(',') + str(blastp_np_2[i,0])
    else:
        blastp_np_2_filtered = np.append(blastp_np_2_filtered, [blastp_np_2[i,]], axis=0)


#convert np array to pd df
blastp_df_2_filtered = pd.DataFrame(blastp_np_2_filtered, columns=["qseqid2", "sseqid2", "pident2", "length2","e-value2"])

with open(str(cwd + '/' + 'blastp_df_2_filtered.out'), 'w') as out_file:
    blastp_df_2_filtered.to_csv(out_file, sep="\t", header=False, index=False)

# modify pfam out: get one record for one transcript (best scored).
#convert pd df to np array bcoz it is easier to work with (for me ;) )
pfam_np = np.array(pfam_df)
print(pfam_df.head())
pfam_np_filtered=np.empty((0, pfam_np.shape[1]))
print(pfam_np.shape[0])

pfam_np_filtered = np.append(pfam_np_filtered, [pfam_np[0,]], axis=0)

for i in progressbar.progressbar(range(1, pfam_np.shape[0],1)):
    if(pfam_np[i,2] == pfam_np_filtered[-1,2]):
        pfam_np_filtered[-1,1] = str(pfam_np_filtered[-1,1]) + str('-') + str(pfam_np[i,1])
        pfam_np_filtered[-1,0] = str(pfam_np_filtered[-1,0]) + str('-') + str(pfam_np[i,0])
    else:
        pfam_np_filtered = np.append(pfam_np_filtered, [pfam_np[i,]], axis=0)

print(pfam_np_filtered.shape)

#convert np array to pd df
pfam_df_filtered_1 = pd.DataFrame(pfam_np_filtered, columns=["target", "accession1", "query_name"])
pfam_df_filtered = pfam_df_filtered_1[["accession1", "query_name"]]

with open(str(cwd + '/' + 'pfam_df_filtered.out'), 'w') as out_file:
    pfam_df_filtered.to_csv(out_file, sep="\t", header=False, index=False)

#to save memory we should delte df which we will not use
del(blastp_df_1)
del(blastp_df_2)
del(pfam_df)
del(blastp_np_1)
del(blastp_np_2)
del(blastp_np_1_filtered)
del(blastp_np_2_filtered)

#merge these data one by one
#first we will merge blastp with pfam

join1_1 = pd.merge(blastp_df_1_filtered, tmhmm_df, how='left', left_on='qseqid', right_on='target_tm')
#print(join1.head())
join1_1_1 = pd.merge(blastp_df_2_filtered, join1_1, how='left', left_on='qseqid2', right_on='qseqid')
#print(join1.head())
#join1_filter = join1[join1['qseqid']!=np.nan]
join2_1 = pd.merge(join1_1_1, signalp_df, how='left', left_on='qseqid', right_on='target_sp')
#print(join2.head(10))
#join2_filter = join2[join2['qseqid']!=np.nan]
#join1_2 = pd.merge(pfam_df_filtered, tmhmm_df, how='left', left_on='query_name', right_on='target_tm')
#print(join3.head(10))
#join3_filter = join3[join3['qseqid']!=np.nan]
#join2_2 = pd.merge(join1_2, signalp_df, how='left', left_on='query_name', right_on='target_sp')
#print(final.head(10))

join3_1_1 = pd.merge(gff_df, join2_1, how='left', left_on='target_gtf', right_on='qseqid')
join3_2 = pd.merge(gff_df, pfam_df_filtered, how='left', left_on='target_gtf', right_on='query_name')

with open(str(cwd + '/' + 'join1_1'), 'w') as out_file:
    join1_1.to_csv(out_file, sep="\t", header=True, index=False)
with open(str(cwd + '/' + 'join1_1_1'), 'w') as out_file:
    join1_1_1.to_csv(out_file, sep="\t", header=True, index=False)
#with open(str(cwd + '/' + 'join1_2'), 'w') as out_file:
#    join1_2.to_csv(out_file, sep="\t", header=True, index=False)
with open(str(cwd + '/' + 'join2_1'), 'w') as out_file:
    join2_1.to_csv(out_file, sep="\t", header=True, index=False)
#with open(str(cwd + '/' + 'join2_2'), 'w') as out_file:
#    join2_2.to_csv(out_file, sep="\t", header=True, index=False)
with open(str(cwd + '/' + 'join3_1'), 'w') as out_file:
    join3_1_1.to_csv(out_file, sep="\t", header=True, index=False)
with open(str(cwd + '/' + 'join3_2'), 'w') as out_file:
    join3_2.to_csv(out_file, sep="\t", header=True, index=False)

#to save memory we should delte df which we will not use
del(signalp_df)
del(join1_1)
#del(join1_2)
del(join2_1)
#del(join2_2)
del(tmhmm_df)

###uncomment below if above script fails with memory error
####from here####
#join3_1 = pd.read_table("/media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/annotation/join3_1", delimiter="\t", comment="#", names=["seq_name","source","feature","start","end","score","strand","frame","attribute","target_gtf","qseqid2","sseqid2","pident2","length2","e-value2","target_tm","Topology","target_sp","SignalP"], header=True, low_memory=False)
#join3_2 = pd.read_table("/media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/annotation/join3_2", delimiter="\t", comment="#", names=["seq_name","source","feature","start","end","score","strand","frame","attribute","target_gtf","accession1","query_name"], header=True, low_memory=False)

#join3_1_1 = pd.read_table("/media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/annotation/join3_1", delimiter="\t", comment="#", header=0, low_memory=False)
#join3_2 = pd.read_table("/media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/annotation/join3_2", delimiter="\t", comment="#", header=0, low_memory=False)


#join3_1 = pd.read_table("/media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/annotation/complete_uniprot_blastoutput.txt", delimiter="\t", comment="#", names=["target_gtf","sseqid","pident","length","e-value"], header=None, low_memory=False)
#join3_2 = pd.read_table("/media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/annotation/final_annotated_v4.1_sorted_modified_temp.gff3", delimiter="\t", comment="#", names=["seq_name","source","feature","start","end","score","strand","frame","attribute","target_gtf"], header=None, low_memory=False)

###remove extra columns to save memory
join3_1 = join3_1_1[["target_gtf","qseqid2","sseqid2","pident2","length2","e-value2","qseqid","sseqid","pident","length","e-value","target_tm","Topology","target_sp","SignalP"]]

print(join3_1.head())
print(join3_2.head())

#join3_1_filtered = join3_1[["target_gtf","qseqid","sseqid","pident","length","e-value"]]
join3_1_noNaN = join3_1[join3_1['target_gtf']!="-"]
print(join3_1.head())
print(join3_1_noNaN.head())
del(join3_1)
del(join3_1_1)
#####till here#####

#final = pd.merge(join3_2, join3_1_filtered_noNaN, how='left', left_on='target_gtf', right_on='target_gtf')
final = pd.merge(join3_2, join3_1_noNaN, how='left', left_on='target_gtf', right_on='target_gtf')

#with open(str(cwd + '/' + 'annotated.gff'), 'w') as out_file:
#    final.to_csv(out_file, sep="\t", header=False, index=False)

with open(str(cwd + '/' + 'final_annotated_v4.3.gff3'), 'w') as out_file:
    final.to_csv(out_file, sep="\t", header=False, index=False)
