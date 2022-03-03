

import progressbar
import argparse
import os
#import csv
import re
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import seaborn as sns

cwd = os.getcwd()

parser = argparse.ArgumentParser(description='give cuffdiff output gene_exp.diff')
parser.add_argument('-diff','--diff', type=str, metavar='', required=True, help='cuffdiff')
parser.add_argument('-list','--list', type=str, metavar='', required=True, help='sig expressed genes')
args = parser.parse_args()

diff_file = args.diff
diff_df = pd.read_table(diff_file, delimiter="\t", comment="#", names=["test_id","gene_id","gene","locus","sample_1","sample_2","status","value_1","value_2","log2fold_change","test_stat","p_value","q_value","significant"], low_memory=False)
gene_list = args.list
genes_df = pd.read_table(gene_list, delimiter="\t", comment="#", names=["gene"], low_memory=False)
#convert to np array
diff_df_1 = np.array(diff_df)

#extract all unique ids
ids = np.unique(diff_df_1[:,0])

ids=ids[100:1000]
genes_df_np = np.array(genes_df)
#create empty array to save all different category genes
'''
steep_down=np.empty((0, diff_df_1.shape[1]))
steep_up=np.empty((0, diff_df_1.shape[1]))
constant=np.empty((0, diff_df_1.shape[1]))
up_after4=np.empty((0, diff_df_1.shape[1]))
up_after8=np.empty((0, diff_df_1.shape[1]))
rest=np.empty((0, diff_df_1.shape[1]))
'''
steep_down=np.empty((0, 11))
steep_up=np.empty((0, 11))
constant=np.empty((0, 11))
up_after4=np.empty((0, 11))
up_after8=np.empty((0, 11))
rest=np.empty((0, 11))

'''
steep_down=[]
steep_up=[]
constant=[]
up_after4=[]
up_after8=[]
rest=[]
'''

#loop for all ids to get all stages one by one and then check for pattern in all those stages

def almost_equal(a,b,x):
    lower=b-((b*x)/100)
    upper=b+((b*x)/100)
    if((a<=upper) and (a>=lower)):
        return True

for i in progressbar.progressbar(genes_df_np[:]):
    all_stages = diff_df_1[diff_df_1[:,1]==i]
    if (float(all_stages[0,7])>float(all_stages[3,7])) and (float(all_stages[3,7])>float(all_stages[5,7])) and (float(all_stages[5,7])>float(all_stages[8,8])):
        all_stages_update = np.append(i, all_stages[:,7].T)
        all_stages_update = np.append(all_stages_update,all_stages[8,8])
        steep_down=np.append(steep_down,all_stages_update.reshape(1,11),axis=0)
        #steep_down=steep_down.append(all_stages_update)
    elif (float(all_stages[0,7]<all_stages[3,7])) and float((all_stages[3,7]<all_stages[5,7])) and float((all_stages[5,7]<all_stages[8,8])):
        all_stages_update = np.append(i, all_stages[:,7].T)
        all_stages_update = np.append(all_stages_update,all_stages[8,8])
        #steep_up=steep_up.append(all_stages_update)
        steep_up=np.append(steep_up,all_stages_update.reshape(1,11),axis=0)
    elif (almost_equal(float(all_stages[0,7]),float(all_stages[3,7]),1)) and (almost_equal(float(all_stages[3,7]),float(all_stages[5,7]),1)) and (almost_equal(float(all_stages[5,7]),float(all_stages[8,8]),1)):
        all_stages_update = np.append(i, all_stages[:,7].T)
        all_stages_update = np.append(all_stages_update,all_stages[8,8])
        #constant=constant.append(all_stages_update)
        constant=np.append(constant,all_stages_update.reshape(1,11),axis=0)
    elif (almost_equal(float(all_stages[0,7]),float(all_stages[3,7]),1)) and (float(all_stages[3,7])<float(all_stages[5,7])) and (float(all_stages[5,7])<=float(all_stages[8,8])):
        all_stages_update = np.append(i, all_stages[:,7].T)
        all_stages_update = np.append(all_stages_update,all_stages[8,8])
        #up_after4=up_after4.append(all_stages_update)
        up_after4=np.append(up_after4,all_stages_update.reshape(1,11),axis=0)
    elif (almost_equal(float(all_stages[0,7]),float(all_stages[3,7]),1)) and (almost_equal(float(all_stages[3,7]),float(all_stages[5,7]),1)) and (float(all_stages[5,7])<float(all_stages[8,8])):
        all_stages_update = np.append(i, all_stages[:,7].T)
        all_stages_update = np.append(all_stages_update,all_stages[8,8])
        #up_after8=up_after8.append(all_stages_update)
        up_after8=np.append(up_after8,all_stages_update.reshape(1,11),axis=0)
    #else:
        #rest=np.append(rest,all_stages,axis=0)

    #print(steep_up.shape)
    #print(steep_down.shape)
    #print(constant.shape)
    #print(up_after4.shape)
    #print(up_after8.shape)

'''
steep_up_np = np.asarray(steep_up)
steep_down_np = np.asarray(steep_down)
constant_np = np.asarray(constant)
up_after4_np = np.asarray(up_after4)
up_after8_np = np.asarray(up_after8)
'''
steep_up_pd = pd.DataFrame(steep_up, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])
#steep_up_pd = steep_up_pd.set_index("gene_id")
steep_down_pd = pd.DataFrame(steep_down, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])
#steep_down_pd = steep_down_pd.set_index("gene_id")
constant_pd = pd.DataFrame(constant, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])
#constant_pd = constant_pd.set_index("gene_id")
up_after4_pd = pd.DataFrame(up_after4, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])
#up_after4_pd = up_after4_pd.set_index("gene_id")
up_after8_pd = pd.DataFrame(up_after8, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])
#up_after8_pd = up_after8_pd.set_index("gene_id")

with open(str(cwd + '/' + 'steep_up.out'), 'w') as out_file:
    steep_up_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '/' + 'steep_down.out'), 'w') as out_file:
    steep_down_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '/' + 'constant.out'), 'w') as out_file:
    constant_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '/' + 'up_after4.out'), 'w') as out_file:
    up_after4_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '/' + 'up_after8.out'), 'w') as out_file:
    up_after8_pd.to_csv(out_file, sep="\t", header=True, index=False)

#### plot ####
# style
plt.style.use('seaborn-darkgrid')
# create a color palette
palette = plt.get_cmap('Set1')
# figure size in inches
figure(figsize=(12,9), dpi=160)

# steep_up_pd
steep_up_pd = steep_up_pd.drop(columns="s7")
steep_up_pd = steep_up_pd.drop(columns="s9")
steep_up_pd_t = steep_up_pd.T
#print (steep_up_pd_t)
steep_up_pd_t.columns = steep_up_pd_t.iloc[0]
steep_up_pd_t = steep_up_pd_t.drop("gene_id")
#print(steep_up_pd_t.columns)
#print(steep_up_pd_t)
num=0
for column in steep_up_pd_t:
    #print(steep_up_pd_t[column])
    num+=1
    plt.plot(['1','2','3','4','5','6','7','8'], steep_up_pd_t[column], marker='', color=palette(num), linewidth=1, alpha=0.9, label=column)
    #plt.
#steep_up_pd_t['gene_id']
image = cwd + '/' + 'steep_up.png'

plt.savefig(image, )
plt.clf()

# steep_down_pd
steep_down_pd = steep_down_pd.drop(columns="s7")
steep_down_pd = steep_down_pd.drop(columns="s9")
steep_down_pd_t = steep_down_pd.T
#print (steep_down_pd_t)
steep_down_pd_t.columns = steep_down_pd_t.iloc[0]
steep_down_pd_t = steep_down_pd_t.drop("gene_id")
#print(steep_down_pd_t.columns)
#print(steep_down_pd_t)
num=0
for column in steep_down_pd_t:
    #print(steep_down_pd_t[column])
    num+=1
    plt.plot(['1','2','3','4','5','6','7','8'], steep_down_pd_t[column], marker='', color=palette(num), linewidth=1, alpha=0.9, label=column)
#steep_down_pd_t['gene_id']
image = cwd + '/' + 'steep_down.png'
plt.savefig(image)
plt.clf()

# up_after4_pd
up_after4_pd = up_after4_pd.drop(columns="s7")
up_after4_pd = up_after4_pd.drop(columns="s9")
up_after4_pd_t = up_after4_pd.T
#print (up_after4_pd_t)
up_after4_pd_t.columns = up_after4_pd_t.iloc[0]
up_after4_pd_t = up_after4_pd_t.drop("gene_id")
#print(up_after4_pd_t.columns)
#print(up_after4_pd_t)
num=0
for column in up_after4_pd_t:
    #print(up_after4_pd_t[column])
    num+=1
    plt.plot(['1','2','3','4','5','6','7','8'], up_after4_pd_t[column], marker='', color=palette(num), linewidth=1, alpha=0.9, label=column)
#up_after4_pd_t['gene_id']
image = cwd + '/' + 'up_after4_pd.png'
plt.savefig(image)
plt.clf()

# up_after8_pd
up_after8_pd = up_after8_pd.drop(columns="s7")
up_after8_pd = up_after8_pd.drop(columns="s9")
up_after8_pd_t = up_after8_pd.T
#print (up_after8_pd_t)
up_after8_pd_t.columns = up_after8_pd_t.iloc[0]
up_after8_pd_t = up_after8_pd_t.drop("gene_id")
#print(up_after8_pd_t.columns)
#print(up_after8_pd_t)
num=0
for column in up_after8_pd_t:
    #print(up_after8_pd_t[column])
    num+=1
    plt.plot(['1','2','3','4','5','6','7','8'], up_after8_pd_t[column], marker='', color=palette(num), linewidth=1, alpha=0.9, label=column)
#up_after8_pd_t['gene_id']
image = cwd + '/' + 'up_after8_pd.png'
plt.savefig(image)
plt.clf()