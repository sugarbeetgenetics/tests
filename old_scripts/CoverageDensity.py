import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv

import sys
cmdargs = str(sys.argv)

input_file = str(sys.argv[1])
interval= int(sys.argv[2])

coverage_data =  np.array(pd.read_csv(input_file, sep='\t'))
i=0
totals = []
bins = []

ids = []

for i in range(1,len(coverage_data[:,1])-1):
    j=0
    total = 0
    for j in range(0,interval):
        if (coverage_data[i,1]==coverage_data[i+1,1]):
            total = coverage_data[i,3]+total
            j=j+1
            i=i+1
        else:
            j=j+1
            i=i+1
            break
        totals = totals.append(total)
        bins = bins.append(j)
        ids = ids.append(coverage_data[i,1])
print (bins,totals,sep='\t',end='\n',file=sys.stdout)

'''
zipped = zip(ids,bins, totals)
with open(str(sys.argv[3]),'w+') as f:
    writer = csv.writer(f,delimiter='\t')
    writer.writerows(zipped)
'''