#usage
#python3 /media/pflanz/Hs1-2/Scripts/combine_annotation_gff_toNewGFF/get_different_pattern_genes_v1.4.py -diff_all /media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/gene_exp_filtered.diff -diff_14 13_gene_exp.diff -diff_47 36_gene_exp.diff -diff_78 68_gene_exp.diff -diff_89 89_gene_exp.diff -diff_79 79_gene_exp.diff -list 13_36_69_sigDiff.genes

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


version="v1_4"
cwd = os.getcwd() + "/" + version

parser = argparse.ArgumentParser(description='give cuffdiff output gene_exp.diff')
parser.add_argument('-diff_all','--diff_all', type=str, metavar='', required=True, help='cuffdiff_all')
parser.add_argument('-diff_14','--diff_14', type=str, metavar='', required=True, help='cuffdiff1_4 q1_q3')
parser.add_argument('-diff_47','--diff_47', type=str, metavar='', required=True, help='cuffdiff4_7 q3_q6')
parser.add_argument('-diff_78','--diff_78', type=str, metavar='', required=True, help='cuffdiff7_8 q6_q8')
parser.add_argument('-diff_89','--diff_89', type=str, metavar='', required=True, help='cuffdiff8_9 q8_q10')
parser.add_argument('-diff_79','--diff_79', type=str, metavar='', required=True, help='cuffdiff7_9 q6_q10')
parser.add_argument('-list','--list', type=str, metavar='', required=True, help='sig expressed genes')
#parser.add_argument('-int_list','--int_list', type=str, metavar='', required=True, help='genes of interest')

args = parser.parse_args()

diff_all_file = args.diff_all
diff_all_df = pd.read_table(diff_all_file, delimiter="\t", comment="#", names=["test_id","gene_id","gene","locus","sample_1","sample_2","status","value_1","value_2","log2fold_change","test_stat","p_value","q_value","significant"], low_memory=False)
diff_14_file = args.diff_14
diff_14_df = pd.read_table(diff_14_file, delimiter="\t", comment="#", names=["test_id","gene_id","gene","locus","sample_1","sample_2","status","value_1","value_2","log2fold_change","test_stat","p_value","q_value","significant"], low_memory=False)
diff_47_file = args.diff_47
diff_47_df = pd.read_table(diff_47_file, delimiter="\t", comment="#", names=["test_id","gene_id","gene","locus","sample_1","sample_2","status","value_1","value_2","log2fold_change","test_stat","p_value","q_value","significant"], low_memory=False)
diff_78_file = args.diff_78
diff_78_df = pd.read_table(diff_78_file, delimiter="\t", comment="#", names=["test_id","gene_id","gene","locus","sample_1","sample_2","status","value_1","value_2","log2fold_change","test_stat","p_value","q_value","significant"], low_memory=False)
diff_89_file = args.diff_89
diff_89_df = pd.read_table(diff_89_file, delimiter="\t", comment="#", names=["test_id","gene_id","gene","locus","sample_1","sample_2","status","value_1","value_2","log2fold_change","test_stat","p_value","q_value","significant"], low_memory=False)
diff_79_file = args.diff_79
diff_79_df = pd.read_table(diff_79_file, delimiter="\t", comment="#", names=["test_id","gene_id","gene","locus","sample_1","sample_2","status","value_1","value_2","log2fold_change","test_stat","p_value","q_value","significant"], low_memory=False)

gene_list = args.list
genes_df = pd.read_table(gene_list, delimiter="\t", comment="#", names=["gene"], low_memory=False)
#int_gene_list = args.int_list
#int_genes = pd.read_table(int_gene_list, delimiter="\t", comment="#", names=["gene_name", "id"], low_memory=False)

#convert to np array
diff_all_df_1 = np.array(diff_all_df)
diff_14_df_1 = np.array(diff_14_df)
diff_47_df_1 = np.array(diff_47_df)
diff_79_df_1 = np.array(diff_79_df)
diff_78_df_1 = np.array(diff_78_df)
diff_89_df_1 = np.array(diff_89_df)
genes_df_np = np.array(genes_df)
#int_genes_np = np.array(int_genes)
#create empty array to save all different category genes

steep_down=np.empty((0, 11))
steep_down_14_only=np.empty((0, 11))
steep_down_14=np.empty((0, 11))
steep_down_47=np.empty((0, 11))
steep_down_79=np.empty((0, 11))
steep_up=np.empty((0, 11))
steep_up_vlTOvh=np.empty((0, 11))
steep_up_vlTOvh_l=np.empty((0, 11))
steep_up_14=np.empty((0, 11))
steep_up_47=np.empty((0, 11))
steep_up_79=np.empty((0, 11))
rest=np.empty((0, 11))

#int_gene_data = np.empty((0,11))

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

for i in tqdm(genes_df_np[:]):
    all_stages = diff_all_df_1[diff_all_df_1[:,1]==i]
    ### down
    if (float(all_stages[0,7])>float(all_stages[2,7])) and (float(all_stages[2,7])>float(all_stages[5,7])) and (float(all_stages[5,7])>float(all_stages[8,8])):
        #steep down 14 47 79
        if (diff_14_df_1[diff_14_df_1[:,1]==i,13]=="yes") and (diff_47_df_1[diff_47_df_1[:,1]==i,13]=="yes") and (diff_79_df_1[diff_79_df_1[:,1]==i,13]=="yes"):
            if (float(all_stages[0,7])<float(all_stages[1,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][0]!="yes")) or (float(all_stages[0,7])>=float(all_stages[1,7])):
                if (float(all_stages[1,7])<float(all_stages[2,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][1]!="yes")) or (float(all_stages[1,7])>=float(all_stages[2,7])):
                    if (float(all_stages[2,7])<float(all_stages[3,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][2]!="yes")) or (float(all_stages[2,7])>=float(all_stages[3,7])):
                        if (float(all_stages[3,7])<float(all_stages[4,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][3]!="yes")) or (float(all_stages[3,7])>=float(all_stages[4,7])):
                            if (float(all_stages[4,7])<float(all_stages[5,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][4]!="yes")) or (float(all_stages[4,7])>=float(all_stages[5,7])):
                                if (float(all_stages[5,7])<float(all_stages[7,7]) and (diff_78_df_1[diff_78_df_1[:,1]==i,13]!="yes")) or (float(all_stages[5,7])>=float(all_stages[7,7])):
                                    if (float(all_stages[7,7])<float(all_stages[8,8]) and (diff_89_df_1[diff_89_df_1[:,1]==i,13]!="yes")) or (float(all_stages[7,7])>=float(all_stages[8,8])):
                                        all_stages_update = np.append(i, all_stages[:,7].T)
                                        all_stages_update = np.append(all_stages_update, all_stages[8,8])
                                        steep_down=np.append(steep_down,all_stages_update.reshape(1,11),axis=0)
        else:
            all_stages_update = np.append(i, all_stages[:,7].T)
            all_stages_update = np.append(all_stages_update, all_stages[8,8])
            rest=np.append(rest,all_stages_update.reshape(1,11),axis=0)

        #steep down 14 only
        if (diff_14_df_1[diff_14_df_1[:,1]==i,13]=="yes") and (diff_47_df_1[diff_47_df_1[:,1]==i,13]!="yes") and (diff_79_df_1[diff_79_df_1[:,1]==i,13]!="yes"):
            if ((float(all_stages[0,7])>=1.5) and (float(all_stages[2,7])>=0 and float(all_stages[2,7])<=1.5) and (float(all_stages[3,7])>=0 and float(all_stages[3,7])<=1.5) and (float(all_stages[4,7])>=0 and float(all_stages[4,7])<=1.5) and (float(all_stages[5,7])>=0 and float(all_stages[5,7])<=1.5)):
                if (float(all_stages[0,7])<float(all_stages[1,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][0]!="yes")) or (float(all_stages[0,7])>=float(all_stages[1,7])):
                    if (float(all_stages[1,7])<float(all_stages[2,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][1]!="yes")) or (float(all_stages[1,7])>=float(all_stages[2,7])):
                        if (float(all_stages[2,7])<float(all_stages[3,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][2]!="yes")) or (float(all_stages[2,7])>=float(all_stages[3,7])):
                            if (float(all_stages[3,7])<float(all_stages[4,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][3]!="yes")) or (float(all_stages[3,7])>=float(all_stages[4,7])):
                                if (float(all_stages[4,7])<float(all_stages[5,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][4]!="yes")) or (float(all_stages[4,7])>=float(all_stages[5,7])):
                                    if (float(all_stages[5,7])<float(all_stages[7,7]) and (diff_78_df_1[diff_78_df_1[:,1]==i,13]!="yes")) or (float(all_stages[5,7])>=float(all_stages[7,7])):
                                        if (float(all_stages[7,7])<float(all_stages[8,8]) and (diff_89_df_1[diff_89_df_1[:,1]==i,13]!="yes")) or (float(all_stages[7,7])>=float(all_stages[8,8])):
                                            all_stages_update = np.append(i, all_stages[:,7].T)
                                            all_stages_update = np.append(all_stages_update, all_stages[8,8])
                                            steep_down_14_only=np.append(steep_down_14_only,all_stages_update.reshape(1,11),axis=0)

        #steep down 14 
        if (diff_14_df_1[diff_14_df_1[:,1]==i,13]=="yes") and (diff_47_df_1[diff_47_df_1[:,1]==i,13]!="yes") and (diff_79_df_1[diff_79_df_1[:,1]==i,13]!="yes"):
            if (float(all_stages[0,7])<float(all_stages[1,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][0]!="yes")) or (float(all_stages[0,7])>=float(all_stages[1,7])):
                if (float(all_stages[1,7])<float(all_stages[2,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][1]!="yes")) or (float(all_stages[1,7])>=float(all_stages[2,7])):
                    if (float(all_stages[2,7])<float(all_stages[3,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][2]!="yes")) or (float(all_stages[2,7])>=float(all_stages[3,7])):
                        if (float(all_stages[3,7])<float(all_stages[4,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][3]!="yes")) or (float(all_stages[3,7])>=float(all_stages[4,7])):
                            if (float(all_stages[4,7])<float(all_stages[5,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][4]!="yes")) or (float(all_stages[4,7])>=float(all_stages[5,7])):
                                if (float(all_stages[5,7])<float(all_stages[7,7]) and (diff_78_df_1[diff_78_df_1[:,1]==i,13]!="yes")) or (float(all_stages[5,7])>=float(all_stages[7,7])):
                                    if (float(all_stages[7,7])<float(all_stages[8,8]) and (diff_89_df_1[diff_89_df_1[:,1]==i,13]!="yes")) or (float(all_stages[7,7])>=float(all_stages[8,8])):
                                        all_stages_update = np.append(i, all_stages[:,7].T)
                                        all_stages_update = np.append(all_stages_update, all_stages[8,8])
                                        steep_down_14=np.append(steep_down_14,all_stages_update.reshape(1,11),axis=0)
        else:
            all_stages_update = np.append(i, all_stages[:,7].T)
            all_stages_update = np.append(all_stages_update, all_stages[8,8])
            rest=np.append(rest,all_stages_update.reshape(1,11),axis=0)

        #steep down 47 
        if (diff_14_df_1[diff_14_df_1[:,1]==i,13]!="yes") and (diff_47_df_1[diff_47_df_1[:,1]==i,13]=="yes"):
            if (float(all_stages[5,7])<float(all_stages[8,8]) and (diff_79_df_1[diff_79_df_1[:,1]==i,13]!="yes")) or float(all_stages[5,7])>=float(all_stages[8,8]):
                if (diff_all_df_1[diff_all_df_1[:,1]==i,13][0]!="yes"):
                    if (diff_all_df_1[diff_all_df_1[:,1]==i,13][1]!="yes"):
                        if (float(all_stages[2,7])<float(all_stages[3,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][2]!="yes")) or (float(all_stages[2,7])>=float(all_stages[3,7])):
                            if (float(all_stages[3,7])<float(all_stages[4,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][3]!="yes")) or (float(all_stages[3,7])>=float(all_stages[4,7])):
                                if (float(all_stages[4,7])<float(all_stages[5,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][4]!="yes")) or (float(all_stages[4,7])>=float(all_stages[5,7])):
                                    if (float(all_stages[5,7])<float(all_stages[7,7]) and (diff_78_df_1[diff_78_df_1[:,1]==i,13]!="yes")) or (float(all_stages[5,7])>=float(all_stages[7,7])):
                                        if (float(all_stages[7,7])<float(all_stages[8,8]) and (diff_89_df_1[diff_89_df_1[:,1]==i,13]!="yes")) or (float(all_stages[7,7])>=float(all_stages[8,8])):
                                            all_stages_update = np.append(i, all_stages[:,7].T)
                                            all_stages_update = np.append(all_stages_update, all_stages[8,8])
                                            steep_down_47=np.append(steep_down_47,all_stages_update.reshape(1,11),axis=0)
        else:
            all_stages_update = np.append(i, all_stages[:,7].T)
            all_stages_update = np.append(all_stages_update, all_stages[8,8])
            rest=np.append(rest,all_stages_update.reshape(1,11),axis=0)

        #steep down 79 
        if (diff_14_df_1[diff_14_df_1[:,1]==i,13]!="yes") and (diff_47_df_1[diff_47_df_1[:,1]==i,13]!="yes"):
            if (diff_79_df_1[diff_79_df_1[:,1]==i,13]=="yes") and float(all_stages[5,7])>float(all_stages[8,8]):
                if (diff_all_df_1[diff_all_df_1[:,1]==i,13][0]!="yes"):
                    if (diff_all_df_1[diff_all_df_1[:,1]==i,13][1]!="yes"):
                        if (diff_all_df_1[diff_all_df_1[:,1]==i,13][2]!="yes"):
                            if (diff_all_df_1[diff_all_df_1[:,1]==i,13][3]!="yes"):
                                if (diff_all_df_1[diff_all_df_1[:,1]==i,13][4]!="yes"):
                                    if (diff_all_df_1[diff_all_df_1[:,1]==i,13][4]!="yes"):
                                        if (float(all_stages[5,7])<float(all_stages[7,7]) and (diff_78_df_1[diff_78_df_1[:,1]==i,13]!="yes")) or (float(all_stages[5,7])>=float(all_stages[7,7])):
                                            if (float(all_stages[7,7])<float(all_stages[8,8]) and (diff_89_df_1[diff_89_df_1[:,1]==i,13]!="yes")) or (float(all_stages[7,7])>=float(all_stages[8,8])):
                                                all_stages_update = np.append(i, all_stages[:,7].T)
                                                all_stages_update = np.append(all_stages_update, all_stages[8,8])
                                                steep_down_79=np.append(steep_down_79,all_stages_update.reshape(1,11),axis=0)
        else:
            all_stages_update = np.append(i, all_stages[:,7].T)
            all_stages_update = np.append(all_stages_update, all_stages[8,8])
            rest=np.append(rest,all_stages_update.reshape(1,11),axis=0)

    ### up
    elif ((float(all_stages[0,7])<float(all_stages[2,7])) and (float(all_stages[2,7])<float(all_stages[5,7])) and (float(all_stages[5,7])<float(all_stages[8,8]))):
        #steep up 14 47 79
        if (diff_14_df_1[diff_14_df_1[:,1]==i,13]=="yes") and (diff_47_df_1[diff_47_df_1[:,1]==i,13]=="yes") and (diff_79_df_1[diff_79_df_1[:,1]==i,13]=="yes"):
            if (float(all_stages[0,7])>float(all_stages[1,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][0]!="yes")) or (float(all_stages[0,7])<=float(all_stages[1,7])):
                if (float(all_stages[1,7])>float(all_stages[2,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][1]!="yes")) or (float(all_stages[1,7])<=float(all_stages[2,7])):
                    if (float(all_stages[2,7])>float(all_stages[3,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][2]!="yes")) or (float(all_stages[2,7])<=float(all_stages[3,7])):
                        if (float(all_stages[3,7])>float(all_stages[4,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][3]!="yes")) or (float(all_stages[3,7])<=float(all_stages[4,7])):
                            if (float(all_stages[4,7])>float(all_stages[5,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][4]!="yes")) or (float(all_stages[4,7])<=float(all_stages[5,7])):
                                if (float(all_stages[5,7])>float(all_stages[7,7]) and (diff_78_df_1[diff_78_df_1[:,1]==i,13]!="yes")) or (float(all_stages[5,7])<=float(all_stages[7,7])):
                                    if (float(all_stages[7,7])>float(all_stages[8,8]) and (diff_89_df_1[diff_89_df_1[:,1]==i,13]!="yes")) or (float(all_stages[7,7])<=float(all_stages[8,8])):
                                        all_stages_update = np.append(i, all_stages[:,7].T)
                                        all_stages_update = np.append(all_stages_update, all_stages[8,8])
                                        steep_up=np.append(steep_up,all_stages_update.reshape(1,11),axis=0)
        else:
            all_stages_update = np.append(i, all_stages[:,7].T)
            all_stages_update = np.append(all_stages_update, all_stages[8,8])
            rest=np.append(rest,all_stages_update.reshape(1,11),axis=0)
        #steep up very low to very high
        if (float(all_stages[8,8])>=10*(float(all_stages[0,7]))):
            if (diff_14_df_1[diff_14_df_1[:,1]==i,13]=="yes") and (diff_47_df_1[diff_47_df_1[:,1]==i,13]=="yes"):
                if (float(all_stages[0,7])>float(all_stages[1,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][0]!="yes")) or (float(all_stages[0,7])<=float(all_stages[1,7])):
                    if (float(all_stages[1,7])>float(all_stages[2,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][1]!="yes")) or (float(all_stages[1,7])<=float(all_stages[2,7])):
                        if (float(all_stages[2,7])>float(all_stages[3,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][2]!="yes")) or (float(all_stages[2,7])<=float(all_stages[3,7])):
                            if (float(all_stages[3,7])>float(all_stages[4,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][3]!="yes")) or (float(all_stages[3,7])<=float(all_stages[4,7])):
                                if (float(all_stages[4,7])>float(all_stages[5,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][4]!="yes")) or (float(all_stages[4,7])<=float(all_stages[5,7])):
                                    if (float(all_stages[5,7])>float(all_stages[7,7]) and (diff_78_df_1[diff_78_df_1[:,1]==i,13]!="yes")) or (float(all_stages[5,7])<=float(all_stages[7,7])):
                                        if (float(all_stages[7,7])>float(all_stages[8,8]) and (diff_89_df_1[diff_89_df_1[:,1]==i,13]!="yes")) or (float(all_stages[7,7])<=float(all_stages[8,8])):
                                            all_stages_update = np.append(i, all_stages[:,7].T)
                                            all_stages_update = np.append(all_stages_update, all_stages[8,8])
                                            steep_up_vlTOvh=np.append(steep_up_vlTOvh,all_stages_update.reshape(1,11),axis=0)
        else:
            all_stages_update = np.append(i, all_stages[:,7].T)
            all_stages_update = np.append(all_stages_update, all_stages[8,8])
            rest=np.append(rest,all_stages_update.reshape(1,11),axis=0)

        #steep up very low to very high but later
        if ((float(all_stages[1,7])>=0 and float(all_stages[1,7])<=1.5) and (float(all_stages[2,7])>=0 and float(all_stages[2,7])<=1.5) and (float(all_stages[3,7])>=0 and float(all_stages[3,7])<=1.5) and (float(all_stages[4,7])>=0 and float(all_stages[4,7])<=1.5) and (float(all_stages[5,7])>=0 and float(all_stages[5,7])<=1.5)):
            if (diff_79_df_1[diff_79_df_1[:,1]==i,13]=="yes") and (diff_14_df_1[diff_14_df_1[:,1]==i,13]!="yes") and (diff_47_df_1[diff_47_df_1[:,1]==i,13]!="yes"):
                if (float(all_stages[0,7])>float(all_stages[1,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][0]!="yes")) or (float(all_stages[0,7])<=float(all_stages[1,7])):
                    if (float(all_stages[1,7])>float(all_stages[2,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][1]!="yes")) or (float(all_stages[1,7])<=float(all_stages[2,7])):
                        if (float(all_stages[2,7])>float(all_stages[3,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][2]!="yes")) or (float(all_stages[2,7])<=float(all_stages[3,7])):
                            if (float(all_stages[3,7])>float(all_stages[4,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][3]!="yes")) or (float(all_stages[3,7])<=float(all_stages[4,7])):
                                if (float(all_stages[4,7])>float(all_stages[5,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][4]!="yes")) or (float(all_stages[4,7])<=float(all_stages[5,7])):
                                    if (float(all_stages[5,7])>float(all_stages[7,7]) and (diff_78_df_1[diff_78_df_1[:,1]==i,13]!="yes")) or (float(all_stages[5,7])<=float(all_stages[7,7])):
                                        if (float(all_stages[7,7])>float(all_stages[8,8]) and (diff_89_df_1[diff_89_df_1[:,1]==i,13]!="yes")) or (float(all_stages[7,7])<=float(all_stages[8,8])):
                                            all_stages_update = np.append(i, all_stages[:,7].T)
                                            all_stages_update = np.append(all_stages_update, all_stages[8,8])
                                            steep_up_vlTOvh_l=np.append(steep_up_vlTOvh_l,all_stages_update.reshape(1,11),axis=0)
        else:
            all_stages_update = np.append(i, all_stages[:,7].T)
            all_stages_update = np.append(all_stages_update, all_stages[8,8])
            rest=np.append(rest,all_stages_update.reshape(1,11),axis=0)

        #steep up 14
        if (diff_14_df_1[diff_14_df_1[:,1]==i,13]=="yes") and (diff_47_df_1[diff_47_df_1[:,1]==i,13]!="yes") and (diff_79_df_1[diff_79_df_1[:,1]==i,13]!="yes"):
            if (float(all_stages[0,7])>float(all_stages[1,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][0]!="yes")) or (float(all_stages[0,7])<float(all_stages[1,7])):
                if (float(all_stages[1,7])>float(all_stages[2,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][1]!="yes")) or (float(all_stages[1,7])<float(all_stages[2,7])):
                    if (float(all_stages[2,7])>float(all_stages[3,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][2]!="yes")) or (float(all_stages[2,7])<float(all_stages[3,7])):
                        if (float(all_stages[3,7])>float(all_stages[4,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][3]!="yes")) or (float(all_stages[3,7])<float(all_stages[4,7])):
                            if (float(all_stages[4,7])>float(all_stages[5,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][4]!="yes")) or (float(all_stages[4,7])<float(all_stages[5,7])):
                                if (float(all_stages[5,7])>float(all_stages[7,7]) and (diff_78_df_1[diff_78_df_1[:,1]==i,13]!="yes")) or (float(all_stages[5,7])<float(all_stages[7,7])):
                                    if (float(all_stages[7,7])>float(all_stages[8,8]) and (diff_89_df_1[diff_89_df_1[:,1]==i,13]!="yes")) or (float(all_stages[7,7])<float(all_stages[8,8])):
                                        all_stages_update = np.append(i, all_stages[:,7].T)
                                        all_stages_update = np.append(all_stages_update, all_stages[8,8])
                                        steep_up_14=np.append(steep_up_14,all_stages_update.reshape(1,11),axis=0)
        else:
            all_stages_update = np.append(i, all_stages[:,7].T)
            all_stages_update = np.append(all_stages_update, all_stages[8,8])
            rest=np.append(rest,all_stages_update.reshape(1,11),axis=0)

        #steep up 47
        if (diff_14_df_1[diff_14_df_1[:,1]==i,13]!="yes") and (diff_47_df_1[diff_47_df_1[:,1]==i,13]=="yes"):
            if (diff_all_df_1[diff_all_df_1[:,1]==i,13][0]!="yes"):
                if (diff_all_df_1[diff_all_df_1[:,1]==i,13][1]!="yes"):
                    if (float(all_stages[2,7])>float(all_stages[3,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][2]!="yes")) or (float(all_stages[2,7])<float(all_stages[3,7])):
                        if (float(all_stages[3,7])>float(all_stages[4,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][3]!="yes")) or (float(all_stages[3,7])<float(all_stages[4,7])):
                            if (float(all_stages[4,7])>float(all_stages[5,7]) and (diff_all_df_1[diff_all_df_1[:,1]==i,13][4]!="yes")) or (float(all_stages[4,7])<float(all_stages[5,7])):
                                if (float(all_stages[5,7])>float(all_stages[7,7]) and (diff_78_df_1[diff_78_df_1[:,1]==i,13]!="yes")) or (float(all_stages[5,7])<float(all_stages[7,7])):
                                    if (float(all_stages[7,7])>float(all_stages[8,8]) and (diff_89_df_1[diff_89_df_1[:,1]==i,13]!="yes")) or (float(all_stages[7,7])<float(all_stages[8,8])):
                                        all_stages_update = np.append(i, all_stages[:,7].T)
                                        all_stages_update = np.append(all_stages_update, all_stages[8,8])
                                        steep_up_47=np.append(steep_up_47,all_stages_update.reshape(1,11),axis=0)
        else:
            all_stages_update = np.append(i, all_stages[:,7].T)
            all_stages_update = np.append(all_stages_update, all_stages[8,8])
            rest=np.append(rest,all_stages_update.reshape(1,11),axis=0)

        #steep up 79
        if (diff_14_df_1[diff_14_df_1[:,1]==i,13]!="yes") and (diff_47_df_1[diff_47_df_1[:,1]==i,13]!="yes") and (diff_79_df_1[diff_79_df_1[:,1]==i,13]=="yes"):
            if (diff_all_df_1[diff_all_df_1[:,1]==i,13][0]!="yes"):
                if (diff_all_df_1[diff_all_df_1[:,1]==i,13][1]!="yes"):
                    if (diff_all_df_1[diff_all_df_1[:,1]==i,13][2]!="yes"):
                        if (diff_all_df_1[diff_all_df_1[:,1]==i,13][3]!="yes"):
                            if (diff_all_df_1[diff_all_df_1[:,1]==i,13][4]!="yes"):
                                if (float(all_stages[5,7])>float(all_stages[7,7]) and (diff_78_df_1[diff_78_df_1[:,1]==i,13]!="yes")) or (float(all_stages[5,7])<float(all_stages[7,7])):
                                    if (float(all_stages[7,7])>float(all_stages[8,8]) and (diff_89_df_1[diff_89_df_1[:,1]==i,13]!="yes")) or (float(all_stages[7,7])<float(all_stages[8,8])):
                                        all_stages_update = np.append(i, all_stages[:,7].T)
                                        all_stages_update = np.append(all_stages_update, all_stages[8,8])
                                        steep_up_79=np.append(steep_up_79,all_stages_update.reshape(1,11),axis=0)
        else:
            all_stages_update = np.append(i, all_stages[:,7].T)
            all_stages_update = np.append(all_stages_update, all_stages[8,8])
            rest=np.append(rest,all_stages_update.reshape(1,11),axis=0)

del(all_stages_update)
'''
### get data for gnees of interest
for i in tqdm(int_genes_np[:,1]):
    all_stages = diff_all_df_1[diff_all_df_1[:,1]==i]
    all_stages_update = np.append(i, all_stages[:,7].T)
    all_stages_update = np.append(all_stages_update, all_stages[8,8])
    int_gene_data=np.append(int_gene_data,all_stages_update.reshape(1,11),axis=0)
'''
steep_up_pd = pd.DataFrame(steep_up, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])
steep_up_vlTOvh_pd = pd.DataFrame(steep_up_vlTOvh, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])
steep_up_vlTOvh_l_pd = pd.DataFrame(steep_up_vlTOvh_l, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])
steep_up_14_pd = pd.DataFrame(steep_up_14, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])
steep_up_47_pd = pd.DataFrame(steep_up_47, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])
steep_up_79_pd = pd.DataFrame(steep_up_79, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])

steep_down_pd = pd.DataFrame(steep_down, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])
steep_down_14_only_pd = pd.DataFrame(steep_down_14_only, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])
steep_down_14_pd = pd.DataFrame(steep_down_14, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])
steep_down_47_pd = pd.DataFrame(steep_down_47, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])
steep_down_79_pd = pd.DataFrame(steep_down_79, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])

rest_pd = pd.DataFrame(rest, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])

#int_gene_data_pd = pd.DataFrame(int_gene_data, columns = ['gene_id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'])

#save all the tables into individual files

with open(str(cwd + '_' + 'steep_up.out'), 'w') as out_file:
    steep_up_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '_' + 'steep_up_vlTOvh.out'), 'w') as out_file:
    steep_up_vlTOvh_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '_' + 'steep_up_vlTOvh_l.out'), 'w') as out_file:
    steep_up_vlTOvh_l_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '_' + 'steep_up_14.out'), 'w') as out_file:
    steep_up_14_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '_' + 'steep_up_47.out'), 'w') as out_file:
    steep_up_47_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '_' + 'steep_up_79.out'), 'w') as out_file:
    steep_up_79_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '_' + 'steep_down.out'), 'w') as out_file:
    steep_down_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '_' + 'steep_down_14.out'), 'w') as out_file:
    steep_down_14_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '_' + 'steep_down_14_only.out'), 'w') as out_file:
    steep_down_14_only_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '_' + 'steep_down_47.out'), 'w') as out_file:
    steep_down_47_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '_' + 'steep_down_79.out'), 'w') as out_file:
    steep_down_79_pd.to_csv(out_file, sep="\t", header=True, index=False)

#with open(str(cwd + '_' + 'int_genes.out'), 'w') as out_file:
#    int_gene_data_pd.to_csv(out_file, sep="\t", header=True, index=False)

with open(str(cwd + '_' + 'rest.out'), 'w') as out_file:
    rest_pd.to_csv(out_file, sep="\t", header=True, index=False)

'''
################
print (int_gene_data)
print(int_gene_data_pd)

int_gene_data_pd = int_gene_data_pd.drop(columns="s7")
int_gene_data_pd = int_gene_data_pd.drop(columns="s9")
int_gene_data_pd_t = int_gene_data_pd.T
print(int_gene_data_pd_t)

int_gene_data_pd_t.columns = int_gene_data_pd_t.iloc[0]
int_gene_data_pd_t = int_gene_data_pd_t.drop("gene_id")
print(int_gene_data_pd_t)
#for column in int_gene_data_pd:
    #log transform all columns
#    int_gene_data_pd_t[column] = pd.to_numeric(int_gene_data_pd_t[column]) + 1
#    int_gene_data_pd_t[column] = np.log2(pd.to_numeric(int_gene_data_pd_t[column]))
################
'''


#define function to plot graph
def plot_graph(in_data,text):
    #### plot ####
    # style
    plt.style.use('seaborn-bright')
    # create a color palette
    #palette = plt.get_cmap('Set1')
    # figure size in inches
    figure(figsize=(12,9), dpi=160)

    ## in order to plot broken y-axis we first need to plot data two times 
    # and them limit axis followed by print either side of both plt together.

    # in_data
    in_data = in_data.drop(columns="s7")
    in_data = in_data.drop(columns="s9")
    in_data_t = in_data.T
    #print (in_data_t)
    in_data_t.columns = in_data_t.iloc[0]
    in_data_t = in_data_t.drop("gene_id")
    #print(in_data_t.columns)
    #print(in_data_t)
    fpkm_break = 20
    max_fpkm = fpkm_break
    for column in in_data_t:
        #log transform all columns
        in_data_t[column] = pd.to_numeric(in_data_t[column]) + 1
        in_data_t[column] = np.log2(pd.to_numeric(in_data_t[column]))
        if(in_data_t[column].max()>max_fpkm):
            max_fpkm=in_data_t[column].max()
        if(max_fpkm>fpkm_break):
            #create two subplots
            f,(plt1, plt2) = plt.subplots(2,1,sharex=True, facecolor='w',gridspec_kw={'height_ratios': [1, 2]}, figsize=(12,9))
            f.tight_layout()
    
    num=0
    for column in in_data_t:
        #print(in_data_t[column])
        num+=1
        #find if there are values more than threshold(= 20)
        
        if(max_fpkm>fpkm_break):
            #create two subplots
            plt1.plot(['1','2','4','5','6','7','8','9'], in_data_t[column], marker='', color='#878787', linewidth=1, alpha=0.9, label='_nolegend_') # add color=palette(num) if colors are needed
            plt2.plot(['1','2','4','5','6','7','8','9'], in_data_t[column], marker='', color='#878787', linewidth=1, alpha=0.9, label='_nolegend_')
        else:
            plt.plot(['1','2','4','5','6','7','8','9'], in_data_t[column], marker='', color='#878787', linewidth=1, alpha=0.9, label='_nolegend_')
    
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
    plt.xlabel("weeks after vernalization")
    plt.ylabel("log2(1+RPKM)")
    plt.savefig(image)
    plt.clf()

plot_graph(steep_down_pd,"steep_down_pd")
plot_graph(steep_down_14_only_pd,"steep_down_14_only_pd")
plot_graph(steep_down_14_pd,"steep_down_14_pd")
plot_graph(steep_down_47_pd,"steep_down_47_pd")
plot_graph(steep_down_79_pd,"steep_down_79_pd")
plot_graph(steep_up_pd,"steep_up_pd")
plot_graph(steep_up_vlTOvh_pd,"steep_up_vlTOvh_pd")
plot_graph(steep_up_vlTOvh_l_pd,"steep_up_vlTOvh_l_pd")
plot_graph(steep_up_14_pd,"steep_up_14_pd")
plot_graph(steep_up_47_pd,"steep_up_47_pd")
plot_graph(steep_up_79_pd,"steep_up_79_pd")