import pandas as pd
import numpy as np

def same_columns(data_frame):
    uni_col=pd.unique(data_frame).shape[0]
    return (int(uni_col))

table="/media/pflanz/Hs1-2/Brassica_Napus/gatk/rnaseq_variants/all_stages.txt"

names=["P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"]

table_data=pd.read_csv(table, delim_whitespace=True, header=None)
nrow=table_data.shape[0]

summary=[]
file="/media/pflanz/Hs1-2/Brassica_Napus/gatk/rnaseq_variants/all_stage_summary.txt"

for row in range(0,nrow,1):
    summary.append(table_data.iloc[row].unique())

summary_df=pd.DataFrame(summary)

summary_df.to_csv(file, index=False, header=False)