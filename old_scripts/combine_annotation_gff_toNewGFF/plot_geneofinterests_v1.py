#import progressbar
import argparse
import os
#import csv
import re
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


version="v1_3"
cwd = os.getcwd() + "/" + version

parser = argparse.ArgumentParser(description='give cuffdiff output gene_exp.diff')
parser.add_argument('-diff_all','--diff_all', type=str, metavar='', required=True, help='cuffdiff_all')
parser.add_argument('-list','--list', type=str, metavar='', required=True, help='genes of interest')
parser.add_argument('-name','--name', type=str, metavar='', required=True, help='outfile name')

args = parser.parse_args()

all_file = args.diff_all
all_df = pd.read_table(all_file, delimiter="\t", comment="#", names=["test_id","gene_id","gene","locus","sample_1","sample_2","status","value_1","value_2","log2fold_change","test_stat","p_value","q_value","significant"], low_memory=False)

gene_list = args.list
genes_df = pd.read_table(gene_list, delimiter="\t", comment="#", names=["gene"], low_memory=False)
genes_df_np = np.array(genes_df)
all_df_np = np.array(all_df)
int_gene_data = np.empty((0,11))

out_name = args.name

### get data for gnees of interest
for i in (genes_df_np[:]):
    all_stages = all_df_np[all_df_np[:,1]==i]
    all_stages_update = np.append(i, all_stages[:,7].T)
    all_stages_update = np.append(all_stages_update, all_stages[8,8])
    int_gene_data=np.append(int_gene_data,all_stages_update.reshape(1,11),axis=0)

int_gene_data_pd = pd.DataFrame(int_gene_data, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])

int_gene_data_pd = int_gene_data_pd.drop(columns="s7")
int_gene_data_pd = int_gene_data_pd.drop(columns="s9")
int_gene_data_pd_t = int_gene_data_pd.T

int_gene_data_pd_t.columns = int_gene_data_pd_t.iloc[0]
int_gene_data_pd_t = int_gene_data_pd_t.drop("gene_id")
#for column in int_gene_data_pd:
    #log transform all columns
#    int_gene_data_pd_t[column] = pd.to_numeric(int_gene_data_pd_t[column]) + 1
#    int_gene_data_pd_t[column] = np.log2(pd.to_numeric(int_gene_data_pd_t[column]))
################

with open(str(cwd + '_' + out_name + ".out"), 'w') as out_file:
    int_gene_data_pd_t.to_csv(out_file, sep="\t", header=True, index=False)

#define function to plot graph
def plot_graph(in_data,text):
    #### plot ####
    # style
    plt.style.use('seaborn-bright')
    # create a color palette
    palette = plt.get_cmap('Set1')
    # figure size in inches
    figure(figsize=(12,9), dpi=160)

    ## in order to plot broken y-axis we first need to plot data two times 
    # and them limit axis followed by print either side of both plt together.
    #print(in_data_t.columns)
    #print(in_data_t)
    fpkm_break = 20
    max_fpkm = fpkm_break
    for column in in_data:
        #log transform all columns
        in_data[column] = pd.to_numeric(in_data[column]) + 1
        in_data[column] = np.log2(pd.to_numeric(in_data[column]))
        if(in_data[column].max()>max_fpkm):
            max_fpkm=in_data[column].max()
        if(max_fpkm>fpkm_break):
            #create two subplots
            f,(plt1, plt2) = plt.subplots(2,1,sharex=True, facecolor='w',gridspec_kw={'height_ratios': [1, 2]}, figsize=(12,9))
            f.tight_layout()
    
    num=0
    for column in in_data:
        #print(in_data_t[column])
        num+=1
        #find if there are values more than threshold(= 20)
        plt.plot(['1','2','4','5','6','7','8','9'], in_data[column], marker='', color=palette(num), linewidth=1, alpha=0.9, label=column)

    image = cwd + '_' + text + ".png"
    plt.legend(loc='upper left') #uncoment if you want to show legends
    plt.title(text)
    plt.xlabel("weeks after vernalization")
    plt.ylabel("log2(1+RPKM)")
    plt.savefig(image)
    plt.clf()

plot_graph(int_gene_data_pd_t,out_name)
