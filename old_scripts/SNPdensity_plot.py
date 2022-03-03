import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import sys
cmdargs = str(sys.argv)

snpden = str(sys.argv[1])
vcf_file= str(sys.argv[1])
coverage = str(sys.argv[1])

def SNPdensity(snpden):
    density_data =  np.array(pd.read_csv(snpden, sep='\t'))
    plt.bar(density_data[:,1],density_data[:,2],color='Red',align='edge',width=1)
    plt.xticks(np.arange(min(density_data[:,1]), max(density_data[:,1]), 10000))
    plt.xticks(rotation=90)
    plt.show()

def SNPquality(vcf_file):
    vcf_data = np.array(pd.read_csv(vcf_file, sep='\t', comment='#'))
    plt.bar(vcf_data[:,0], vcf_data[:,6], color='Blue')
    plt.xticks(np.arange(min(vcf_data[:,1]), max(vcf_data[:,1])+1000, 10000))
    plt.show()

def Coverage(coverage):
    density_data =  np.array(pd.read_csv(coverage, sep='\t'))
    plt.bar(density_data[:,1],density_data[:,2],color='Red',align='edge',width=1)
    #plt.xticks(np.arange(min(density_data[:,0]), max(density_data[:,1]), 10000))
    plt.xticks(rotation=90)
    plt.show()

if (str(sys.argv[2])=="snpden"):
    SNPdensity(snpden)
elif (str(sys.argv[2])=="quality"):
    SNPquality(vcf_file)
elif (str(sys.argv[2])=="coverage"):
    Coverage(coverage)