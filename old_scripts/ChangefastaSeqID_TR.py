from os import path
import sys, os, glob
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

parser = argparse.ArgumentParser(description='Change fasta seq id.')
parser.add_argument('-in_file', '--list', type=str, metavar='', required=True, help='squence file, fasta format')
parser.add_argument('-name_list', '--file', type=str, metavar='', required=True, help='squence file, fasta format')
parser.add_argument('-sep', '--sep', type=str, metavar='', required=False, help='name seperator')
parser.add_argument('-sname', '--sample', type=int, metavar='', required=True, help='sample name')


args = parser.parse_args()
fasta_file = args.file

list_file = args.list

sep = args.sep

sample_name = args.sample

cwd = os.getcwd()

list_names = pd.DataFrame(pd.read_csv(list_file, sep = '\t', header = None))

with open(fasta_file) as fasta:
    for seq in SeqIO.parse(fasta, "fasta"):
        old_ID = str(seq.id)
        #id_Strip = old_ID.split(" ")  #will return a list
        #new_ID = paste(sample_name, id_Strip[0], sep="_")
        
        new_sequence = SeqRecord(Seq(str(seq.seq)), id=new_ID, description="")
        with open(str(cwd+'/'+sample_name+'changed_id.fasta'), 'a') as newfile:
            SeqIO.write(new_sequence, newfile, "fasta")
