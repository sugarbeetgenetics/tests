'''
Author: Avneesh Kumar
At Plant Breeding Institute, CAU Kiel

Change headers for all the sequences in a file by adding file name in the sequence header. Provide directory with sequence fasta files.")
'''

from os import path
import sys, os, glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from progressbar import ProgressBar

pbar = ProgressBar()

args = sys.argv

file1 = "/media/planz/Hs1-2/Beet_translocation/UnorderedDataSet/P_procumbens/P_procumbens-1.0.soap.1.fa"
new_file = "/media/planz/Hs1-2/Beet_translocation/UnorderedDataSet/P_procumbens/P_procumbens_1_soap_newNames.fa"

'''
with open(output_file, 'w') as w_file:
    for fil in pbar(fils):
        with open(fil, 'r'):
            seq_records = SeqIO.parse(fil,'fasta')
            SeqIO.write(seq_records, w_file,'fasta')
'''
with open(new_file, 'w') as w_file:
	for seq in SeqIO.parse(file1, "fasta"):
		old_ID = str(seq.id)
		new_ID = "%s%s%s" %(str("Proc_seq1_soap"),str("_"),str(old_ID))
		new_sequence = SeqRecord(Seq(str(seq.seq)), id = new_ID, description = "")
		SeqIO.write(new_sequence, w_file,"fasta")


