#! /usr/bin/python
'''
To extract sequence from SPAdes fasta file by giving their ID in a list

Author: Avneesh Kumar
'''


from progressbar import ProgressBar
pbar=ProgressBar()
import argparse
from os import path
import sys, os, glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

parser = argparse.ArgumentParser(description='This script can filter fasta file based on certain criterias. --written by Avneesh Kumar')
parser.add_argument('-f','--fasta_file',help='fasta file', required=True)
parser.add_argument('-l','--length',help='minimum length of fasta sequence. Default 1000', type=int, default=1000)
parser.add_argument('-N','--percent_N',help='Percent Ns. Default 100', required=False, type=int, default=100)
parser.add_argument('-c','--cov',help='Minimum kmer coverage. Default 10', required=False, type=int, default=10)
parser.add_argument('-o','--output',help='output suffix', required=True)

args = parser.parse_args()
###correct input file and output path
fasta_file = args.fasta_file
print ('Filtering fasta started\n')

i=0

output_path = os.path.join(os.getcwd(),str(args.output)+".fasta")

open(output_path, 'w').close() ##to clean the file if already exist.

for record in pbar(SeqIO.parse(fasta_file, "fasta")):
    if args.fasta_file!=output_path:
        with open(output_path, 'a') as w_file:
            record_id = str(record.id)
            cov = float(record_id.split("_")[5])
            n = record.seq.count("N")
            percent_n = n/len(record.seq)
            if (len(record.seq)>args.length) and percent_n<args.percent_N and cov>=args.cov:
                SeqIO.write(record, w_file,"fasta")
            else:
                pass

print ('Filtering Done!')
