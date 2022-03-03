#!/bin/python3

import argparse
import numpy as np
import os
import pandas as pd

snp_file = pd.read_csv("/media/pflanz/Avneesh/Sarah_FlapJack/final.sorted.snpFile", sep="\t", comment="#", low_memory=False, header=None)

new_rows = pd.DataFrame()

for i,id in enumerate(snp_file.iloc[:,1].unique()):
    all_SRRsIDs = snp_file.iloc[:,0].unique()
    all_SRRsIDs_sort = list(np.sort(all_SRRsIDs))
    subset=snp_file[snp_file.iloc[:,1]==id]
    if i==0:
        new_rows = subset

    extra_id = []
    if subset.shape[0]<25:
        subset_uniq = list(subset.iloc[:,0].unique())
        for item in all_SRRsIDs_sort:
            if item not in subset_uniq:
                extra_id.append(item)
    
    if i!=0:
        new_rows=pd.DataFrame(np.vstack((new_rows, subset)))

    for item in extra_id:
        new_rows=pd.DataFrame(np.vstack((new_rows,np.hstack(([item],subset.iloc[0,1:5],subset.iloc[0,4])))))
    
    #append for the express as one of the accession
    new_rows=pd.DataFrame(np.vstack((new_rows,np.hstack((['Express'],subset.iloc[0,1:5],subset.iloc[0,4])))))

    if(i==4000):
        print("4000 Done")

with open(str('catVCFtoSNPchip_test.txt'), 'w') as out_file:
    new_rows.to_csv(out_file, sep="\t", index=False)
print("new_rows Done")
#new_rows = pd.read_csv("/media/pflanz/Avneesh/catVCFtoSNPchip_test.txt", sep="\t", comment="#", low_memory=False, header=0)

All_IDs = list(new_rows.iloc[:,0].unique())
snpChip = pd.DataFrame(All_IDs)

for i,id in enumerate(new_rows.iloc[:,1].unique()):
    subset=new_rows[new_rows.iloc[:,1]==id]
    snpChip[str(id)] = list(subset.iloc[:,5])
print("snpChip Done")
#header_row = list(new_rows.iloc[:,1].unique()).append(['Accessions'])
#snpChip = pd.DataFrame(np.vstack((header_row,snpChip)))


with open(str('SNPchip_test.txt'), 'w') as out_file:
    snpChip.to_csv(out_file, sep="\t", index=False)

snpChip_HotConvert = pd.DataFrame().reindex_like(snpChip) #create empty df with same structure as snpChip so that we can replace the values
# convert snpChip to 0/1 based on allele (0: Ref, 1: Alt)
for col,id in enumerate(snpChip.iloc[0,:].unique()):
    for row in range(0,27):
        if snpChip.iloc[row,col]==snpChip[26,col]:
            snpChip_HotConvert.iloc[row,col] = 0
        else:
            snpChip_HotConvert.iloc[row,col] = 1
        row+=1
print("snpChip_HotConvert Done")
  
with open(str('SNPchip_HotConvert.txt'), 'w') as out_file:
    snpChip_HotConvert.to_csv(out_file, sep="\t", index=False)
      