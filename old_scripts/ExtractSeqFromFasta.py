'''
To extract sequence from any fasta file by giving their ID in a list
Author: Avneesh Kumar
'''

import argparse
#from progressbar import ProgressBar
#pbar=ProgressBar()
import sys
import sys, os, glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
parser = argparse.ArgumentParser(description='This script can extract sequences from fasta file. --Avneesh Kumar')
parser.add_argument('-f','--fasta_file',help='complete path for fasta file', required=True, type=argparse.FileType('r'))
parser.add_argument('-id','--seq_id',help='Sequence IDs separated by comma(no space is allowed)', required=True)
parser.add_argument('-out','--out_file',help='out file name', required=True)

args = parser.parse_args()

file_name = str(args.out_file)

'''
input_fasta = raw_input('Subject fasta full path:')
ID_list=raw_input('Query ID:')
'''
print ('Reading Subject fasta.... \n')
record_dict = SeqIO.to_dict(SeqIO.parse(args.fasta_file,'fasta'))
print ('Extracting.... \n')
'''
ID_file = open(ID_list,"r")
IDs = ID_file.readlines()
'''
IDs = (str(args.seq_id)).split(",")
cwd = os.getcwd()

i=0
for i in range(len(IDs)):
	outfile = "%s%s%s%s" %(str(cwd),str('/'),file_name,str('.fasta'))
	with open(outfile, 'a') as w_file:
		if (IDs[i] in record_dict.keys()): #check if sequence is present in the fasta file or not.
			SeqIO.write((record_dict[IDs[i]]), w_file,'fasta')# use any record ID
		else:
			print("%s is not present" %(IDs[i]) )
		i +=1
print ('Extraction Done!')
