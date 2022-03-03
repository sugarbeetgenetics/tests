#python3 /media/pflanz/Hs1-2/Scripts/combine_annotation_gff_toNewGFF/get_different_pattern_genes_v1.8.py -diff_all /media/pflanz/BackUps/Brassicanapus/gatk/second_Attempt/gene_exp_filtered.diff --fpkm_all /media/pflanz/BackUps/Brassicanapus/gatk/second_Attempt/gene_exp_analysis/cuffdiff/all_stages_restricted/genes.fpkm_tracking -diff_14 /media/pflanz/BackUps/Brassicanapus/gatk/second_Attempt/gene_exp_analysis/13_gene_exp.diff -diff_47 /media/pflanz/BackUps/Brassicanapus/gatk/second_Attempt/gene_exp_analysis/36_gene_exp.diff -diff_78 /media/pflanz/BackUps/Brassicanapus/gatk/second_Attempt/gene_exp_analysis/78_gene_exp.diff -diff_89 /media/pflanz/BackUps/Brassicanapus/gatk/second_Attempt/gene_exp_analysis/89_gene_exp.diff -diff_79 /media/pflanz/BackUps/Brassicanapus/gatk/second_Attempt/gene_exp_analysis/68_gene_exp.diff

#usage
#python3 /media/pflanz/Hs1-2/Scripts/combine_annotation_gff_toNewGFF/get_different_pattern_genes_v1.6.py 
# --diff_all /media/pflanz/BackUps/Brassicanapus/gatk/second_Attempt/gene_exp_filtered.diff 
# --fpkm_all /media/pflanz/BackUps/Brassicanapus/gatk/second_Attempt/gene_exp_analysis/cuffdiff/all_stages_restricted/genes.fpkm_tracking 
# --diff_14 /media/pflanz/BackUps/Brassicanapus/gatk/second_Attempt/gene_exp_analysis/13_gene_exp.diff 
# --diff_47 /media/pflanz/BackUps/Brassicanapus/gatk/second_Attempt/gene_exp_analysis/36_gene_exp.diff 
# --diff_78 /media/pflanz/BackUps/Brassicanapus/gatk/second_Attempt/gene_exp_analysis/78_gene_exp.diff 
# --diff_89 /media/pflanz/BackUps/Brassicanapus/gatk/second_Attempt/gene_exp_analysis/89_gene_exp.diff 
# --diff_79 /media/pflanz/BackUps/Brassicanapus/gatk/second_Attempt/gene_exp_analysis/68_gene_exp.diff

from tqdm import tqdm
import argparse
import os
import re
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


version="v1_7"
cwd = os.getcwd() + "/" + version

parser = argparse.ArgumentParser(description='give cuffdiff output gene_exp.diff')
parser.add_argument('-diff_all','--diff_all', type=str, metavar='', required=True, help='cuffdiff_all')
parser.add_argument('-fpkm_all','--fpkm_all', type=str, metavar='', required=True, help='fpkm_tracking all')
parser.add_argument('-diff_14','--diff_14', type=str, metavar='', required=True, help='cuffdiff1_4 q1_q3')
parser.add_argument('-diff_47','--diff_47', type=str, metavar='', required=True, help='cuffdiff4_7 q3_q6')
parser.add_argument('-diff_78','--diff_78', type=str, metavar='', required=True, help='cuffdiff7_8 q6_q8')
parser.add_argument('-diff_89','--diff_89', type=str, metavar='', required=True, help='cuffdiff8_9 q8_q10')
parser.add_argument('-diff_79','--diff_79', type=str, metavar='', required=True, help='cuffdiff7_9 q6_q10')

args = parser.parse_args()

fpkm_all_file = args.fpkm_all
fpkm_all_df = pd.read_csv(fpkm_all_file, delimiter="\t", comment="#", header=0)
diff_all_file = args.diff_all
diff_all_df = pd.read_csv(diff_all_file, delimiter="\t", comment="#", names=["test_id","gene_id","gene","locus","sample_1","sample_2","status","value_1","value_2","log2fold_change","test_stat","p_value","q_value","significant"], low_memory=False)
diff_14_file = args.diff_14
diff_14_df = pd.read_csv(diff_14_file, delimiter="\t", comment="#", names=["test_id","gene_id","gene","locus","sample_1","sample_2","status","value_1","value_2","log2fold_change","test_stat","p_value","q_value","significant"], low_memory=False)
diff_47_file = args.diff_47
diff_47_df = pd.read_csv(diff_47_file, delimiter="\t", comment="#", names=["test_id","gene_id","gene","locus","sample_1","sample_2","status","value_1","value_2","log2fold_change","test_stat","p_value","q_value","significant"], low_memory=False)
diff_78_file = args.diff_78
diff_78_df = pd.read_csv(diff_78_file, delimiter="\t", comment="#", names=["test_id","gene_id","gene","locus","sample_1","sample_2","status","value_1","value_2","log2fold_change","test_stat","p_value","q_value","significant"], low_memory=False)
diff_89_file = args.diff_89
diff_89_df = pd.read_csv(diff_89_file, delimiter="\t", comment="#", names=["test_id","gene_id","gene","locus","sample_1","sample_2","status","value_1","value_2","log2fold_change","test_stat","p_value","q_value","significant"], low_memory=False)
diff_79_file = args.diff_79
diff_79_df = pd.read_csv(diff_79_file, delimiter="\t", comment="#", names=["test_id","gene_id","gene","locus","sample_1","sample_2","status","value_1","value_2","log2fold_change","test_stat","p_value","q_value","significant"], low_memory=False)

#convert to np array
diff_all_df_1 = np.array(diff_all_df)
diff_14_df_1 = np.array(diff_14_df)
diff_47_df_1 = np.array(diff_47_df)
diff_79_df_1 = np.array(diff_79_df)
diff_78_df_1 = np.array(diff_78_df)
diff_89_df_1 = np.array(diff_89_df)

#create empty array to save all different category genes

steep_down_14_only=np.empty((0, 9))
steep_up_vlTOvh=np.empty((0, 9))
steep_up_vlTOvh_l=np.empty((0, 9))

#loop for all ids to get all stages one by one and then check for pattern in all those stages

def almost_equal(a,b,x):
    lower1=b-((b*x)/100)
    upper1=b+((b*x)/100)
    lower2=a-((a*x)/100)
    upper2=a+((a*x)/100)
    if((a<=upper1) and (a>=lower1)) or ((b<=upper2) and (b>=lower2)) or (a==b):
        return True
    else:
        return False

total_expressed = 0
for i in tqdm(list(fpkm_all_df['tracking_id'].unique())):
    #print (i)
    all_stages = fpkm_all_df.loc[fpkm_all_df['tracking_id']==i,['tracking_id', 'q1_FPKM','q2_FPKM','q3_FPKM','q4_FPKM','q5_FPKM','q6_FPKM','q7_FPKM','q8_FPKM']]
    all_stages_list = list(all_stages.iloc[0,:])
    #print(all_stages)
    #print(np.array(all_stages_list).reshape(1,9))
    value_count = 0
    for value in all_stages_list[1:]:
        if float(value) >= 1:
            value_count += 1
    
    if value_count == 0:
        total_expressed += 1
        pass

    ### down
    #steep down 14 only
    # genes that are expressed in the vegetative meristem(q1 to q4 or week1 to week 5), and significantly down regulated at q4 or week 5, and remains same as q4 (no positive significant change)
    if (diff_14_df_1[diff_14_df_1[:,1]==i,13]=="yes"): # significant between q1 q4
        if float(all_stages_list[1])>=2.0 and float(all_stages_list[5])<1.0 and float(all_stages_list[6])<=1.0 and float(all_stages_list[7])<=1.0 and float(all_stages_list[8])<=1.0:# set lower/upper limits 
            if (float(all_stages_list[1])<float(all_stages_list[2]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][0]!="yes")) or (float(all_stages_list[1])>=float(all_stages_list[2])):
                if (float(all_stages_list[2])<float(all_stages_list[2]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][1]!="yes")) or (float(all_stages_list[2])>=float(all_stages_list[3])):
                    if (float(all_stages_list[3])<float(all_stages_list[4]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][2]!="yes")) or (float(all_stages_list[3])>=float(all_stages_list[4])):
                        if (float(all_stages_list[4])<float(all_stages_list[5]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][3]!="yes")) or (float(all_stages_list[4])>=float(all_stages_list[5])):
                            if (float(all_stages_list[5])<float(all_stages_list[6]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][4]!="yes")) or (float(all_stages_list[5])>=float(all_stages_list[6])):
                                steep_down_14_only=np.append(steep_down_14_only,np.array(all_stages_list).reshape(1,9),axis=0)

    #steep up very low to very high
    if float(all_stages_list[1])<=10 and float(all_stages_list[4])>=1 and float(all_stages_list[5])>=1 and float(all_stages_list[6])>=1 and float(all_stages_list[7])>=1 and float(all_stages_list[8])>=1: #if (float(all_stages_list[8])>=5*(float(all_stages_list[1]))) and 
        if (diff_14_df_1[diff_14_df_1[:,1]==i,13]=="yes") and (diff_47_df_1[diff_47_df_1[:,1]==i,13]=="yes") and (diff_79_df_1[diff_79_df_1[:,1]==i,13]=="yes"):
            if (float(all_stages_list[1])>float(all_stages_list[2]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][0]!="yes")) or (float(all_stages_list[1])<=float(all_stages_list[2])):
                if (float(all_stages_list[2])>float(all_stages_list[3]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][1]!="yes")) or (float(all_stages_list[2])<=float(all_stages_list[3])):
                    if (float(all_stages_list[3])>float(all_stages_list[4]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][2]!="yes")) or (float(all_stages_list[3])<=float(all_stages_list[4])):
                        if (float(all_stages_list[4])>float(all_stages_list[5]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][3]!="yes")) or (float(all_stages_list[4])<=float(all_stages_list[5])):
                            if (float(all_stages_list[5])>float(all_stages_list[6]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][4]!="yes")) or (float(all_stages_list[5])<=float(all_stages_list[6])):
                                if (float(all_stages_list[6])>float(all_stages_list[7]) and (diff_78_df_1[diff_78_df_1[:,1]==i,13]!="yes")) or (float(all_stages_list[6])<=float(all_stages_list[7])):
                                    if (float(all_stages_list[7])>float(all_stages_list[8]) and (diff_89_df_1[diff_89_df_1[:,1]==i,13]!="yes")) or (float(all_stages_list[7])<=float(all_stages_list[8])):
                                        steep_up_vlTOvh=np.append(steep_up_vlTOvh,np.array(all_stages_list).reshape(1,9),axis=0)
    
    #steep up very low to very high but later
    if float(all_stages_list[1])<=1.0 and float(all_stages_list[2])<=1.0 and float(all_stages_list[3])<=1.0 and float(all_stages_list[4])<=1.0 and float(all_stages_list[5])<=1.0 and float(all_stages_list[8])>=2:
        if (diff_14_df_1[diff_14_df_1[:,1]==i,13]!="yes") and (diff_47_df_1[diff_47_df_1[:,1]==i,13]!="yes"): #(diff_79_df_1[diff_79_df_1[:,1]==i,13]=="yes") and 
            if (float(all_stages_list[1])>float(all_stages_list[2]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][0]!="yes")) or (float(all_stages_list[1])<=float(all_stages_list[2])):
                if (float(all_stages_list[2])>float(all_stages_list[3]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][1]!="yes")) or (float(all_stages_list[2])<=float(all_stages_list[3])):
                    if (float(all_stages_list[3])>float(all_stages_list[4]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][2]!="yes")) or (float(all_stages_list[3])<=float(all_stages_list[4])):
                        if (float(all_stages_list[4])>float(all_stages_list[5]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][3]!="yes")) or (float(all_stages_list[4])<=float(all_stages_list[5])):
                            if (float(all_stages_list[5])>float(all_stages_list[6]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][4]!="yes")) or (float(all_stages_list[5])<=float(all_stages_list[6])):
                                if (float(all_stages_list[6])>float(all_stages_list[7]) and (diff_78_df_1[diff_78_df_1[:,1]==i,13]!="yes")) or (float(all_stages_list[6])<=float(all_stages_list[7])):
                                    if (float(all_stages_list[7])>float(all_stages_list[8]) and (diff_89_df_1[diff_89_df_1[:,1]==i,13]!="yes")) or (float(all_stages_list[7])<=float(all_stages_list[8])):
                                        steep_up_vlTOvh_l=np.append(steep_up_vlTOvh_l,np.array(all_stages_list).reshape(1,9),axis=0)


steep_up_vlTOvh_pd = pd.DataFrame(steep_up_vlTOvh, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8'])
steep_up_vlTOvh_l_pd = pd.DataFrame(steep_up_vlTOvh_l, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8'])
steep_down_14_only_pd = pd.DataFrame(steep_down_14_only, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8'])

#save all the tables into individual files
with open(str(cwd + '_' + 'steep_up_vlTOvh.out'), 'w') as out_file:
    steep_up_vlTOvh_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '_' + 'steep_up_vlTOvh_l.out'), 'w') as out_file:
    steep_up_vlTOvh_l_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '_' + 'steep_down_14_only.out'), 'w') as out_file:
    steep_down_14_only_pd.to_csv(out_file, sep="\t", header=True, index=False)


#define function to plot graph
def plot_graph(in_data, text):
    #### plot ####
    # style
    #plt.style.use('seaborn-bright')
    # create a color palette
    #palette = plt.get_cmap('Set1')
    # figure size in inches
    figure(figsize=(12,9), dpi=160)

    ## in order to plot broken y-axis we first need to plot data two times 
    # and them limit axis followed by print either side of both plt together.

    # in_data
    in_data_t = in_data.T
    #print (in_data_t)
    #in_data_t.columns = in_data_t.iloc[0]
    #in_data_t = in_data_t.drop("gene_id")
    #print(in_data_t.columns)
    #print(in_data_t)
    fpkm_break = 20
    max_fpkm = fpkm_break

    for column in in_data_t.drop("gene_id"):
        #log transform all columns
        in_data_t.iloc[1:,column] = pd.to_numeric(in_data_t.iloc[1:,column]) + 1
        in_data_t.iloc[1:,column] = np.log2(pd.to_numeric(in_data_t.iloc[1:,column]))
        if(in_data_t.iloc[1:,column].max()>max_fpkm):
            max_fpkm=in_data_t.iloc[1:,column].max()
        if(max_fpkm>fpkm_break):
            #create two subplots
            f,(plt1, plt2) = plt.subplots(2,1,sharex=True, facecolor='w',gridspec_kw={'height_ratios': [1, 2]}, figsize=(12,9))
            f.tight_layout()
    
    num=0
    for column in in_data_t.drop("gene_id"):
        #print(in_data_t[column])
        num+=1
        #find if there are values more than threshold(= 20)
        if in_data_t.iloc[0,column] in ["g22155","g52652","g165514"]:
            colour='#878787'
        else:
            colour='#878787'
        if(max_fpkm>fpkm_break):
            #print("here")
            #create two subplots
            plt1.plot(['1','2','4','5','6','7','8','9'], in_data_t.iloc[1:,column], marker='', color=colour, linewidth=1, alpha=0.9, label='_nolegend_') # add color=palette(num) if colors are needed
            ax = plt.gca()
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(1)
                ax.spines[axis].tick_params(width=1)
        
            plt2.plot(['1','2','4','5','6','7','8','9'], in_data_t.iloc[1:,column], marker='', color=colour, linewidth=1, alpha=0.9, label='_nolegend_')
            ax = plt.gca()
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(1)
                ax.spines[axis].tick_params(width=1)
        else:
            
            plt.plot(['1','2','4','5','6','7','8','9'], in_data_t.iloc[1:,column], marker='', color=colour, linewidth=1, alpha=0.9, label='_nolegend_')
            ax = plt.gca()
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(1)
                ax.tick_params(width=1)
    
    # genes of interest
    #int_gene_data_pd = int_gene_data_pd.drop(columns="s7")
    #int_gene_data_pd = int_gene_data_pd.drop(columns="s9")
    #int_gene_data_pd_t = int_gene_data_pd.T
    
    #int_gene_data_pd_t.columns = int_gene_data_pd_t.iloc[0]
    #int_gene_data_pd_t.columns = int_genes_np[:,0]
    #print (int_gene_data_pd_t)
    #int_gene_data_pd_t = int_gene_data_pd_t.drop('gene_id')
    '''
    for column in int_gene_data_pd_t:
        #log transform all columns
        int_gene_data_pd_t[column] = pd.to_numeric(int_gene_data_pd_t[column]) + 1
        int_gene_data_pd_t[column] = np.log2(pd.to_numeric(int_gene_data_pd_t[column]))
    '''
    #    for column in int_gene_data_pd_t:
    #      plt.plot(['1','2','4','5','6','7','8','9'], int_gene_data_pd_t[column], marker='', linewidth=1.5, alpha=0.9, label=column)

    #in_data_t['gene_id']
    image = cwd + '_' + text + ".png"
    if(max_fpkm>fpkm_break): 
        #limit axis
        plt2.set_ylim(0,fpkm_break)
        plt2.yaxis.set_major_locator(plt.MultipleLocator(1.0))
        plt1.set_ylim(fpkm_break,max_fpkm)
        plt1.yaxis.set_major_locator(plt.MultipleLocator(1.0))
        # hide spines between ax and ax2
        plt2.spines['right'].set_visible(False)
        plt1.spines['left'].set_visible(False)
    #plt.legend(loc='upper left') #uncoment if you want to show legends
    plt.title(text)
    plt.xlabel("Weeks during vernalization", fontsize=14) #, fontweight='bold'
    plt.ylabel("log2(1+RPKM)", fontsize=14) #, fontweight='bold'
    plt.savefig(image)
    plt.clf()


plot_graph(steep_down_14_only_pd,"steep_down_14_only_pd")
plot_graph(steep_up_vlTOvh_pd,"steep_up_vlTOvh_pd")
plot_graph(steep_up_vlTOvh_l_pd,"steep_up_vlTOvh_l_pd")
print(total_expressed)