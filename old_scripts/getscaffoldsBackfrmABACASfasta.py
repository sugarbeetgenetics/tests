#######################################################################
#### If A and B are two datasets for sequence scaffolds,       ########
#### this script will give scaffolds present in A but not in B.########
#######################################################################


from progressbar import ProgressBar
pbar=ProgressBar()

from os import path
import sys, os, glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

input_fasta = raw_input('Subject fasta full path:')
ID_list=raw_input('Query ID list full path:')

print ('Reading Subject fasta.... \n')

