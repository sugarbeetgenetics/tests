import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import progressbar

import sys
cmdargs = str(sys.argv)

stage_name=str(sys.argv[1])

vcf_table="/media/pflanz/Hs1-2/Brassica_Napus/gatk/rnaseq_variants/gatk_hc_express_rna_all_stages.vcf.table"

vcf_table_data=pd.read_csv(vcf_table, sep='\t')

def same_columns(data_frame):
    uni_col=pd.unique(data_frame).shape[0]
    total_col=data_frame.shape[0]
    return (total_col-uni_col+1)

file_name="/media/pflanz/Hs1-2/Brassica_Napus/gatk/rnaseq_variants/"+stage_name

vcf_table_data.head()
nrow=vcf_table_data.shape[0]
PASS_FAIL = pd.DataFrame()

exec("%s=[]" % (stage_name))
print(stage_name)
for row in range(0,nrow,1):
    stage1=stage_name+"-"+str(1)+"Aligned.GT"
    stage2=stage_name+"-"+str(2)+"Aligned.GT"
    stage3=stage_name+"-"+str(3)+"Aligned.GT"
    row_subset=vcf_table_data.iloc[row]
    column_subset=row_subset[[stage1,stage2,stage3]]
    if (same_columns(column_subset)>=2):
        exec("%s.append(\"PASS\")" % (stage_name))
    else:
        exec("%s.append(\"FAIL\")" % (stage_name))
        
    print(stage_name,"\t",row,end='\r')

exec("PASS_FAIL['%s']=%s" % (stage_name,stage_name))

PASS_FAIL.to_csv(file_name, index=False, header=False)