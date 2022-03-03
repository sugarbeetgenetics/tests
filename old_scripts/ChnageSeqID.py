'''
Author: Avneesh Kumar
Change headers for all the sequences in a file by adding file name in the sequence header. Provide directory with sequence fasta files.")
'''

from os import path
import sys, os, glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


#from progressbar import ProgressBar
#pbar=ProgressBar()

args = sys.argv

input = args[1]
files = glob.glob("%s/*.f*" %input)

#get the file

'''
with open(output_file, 'w') as w_file:
    for fil in pbar(fils):
        with open(fil, 'r'):
            seq_records = SeqIO.parse(fil,'fasta')
            SeqIO.write(seq_records, w_file,'fasta')
'''

for file in files:
    
    file_name = str(file).split("/")[-1]
    file_name = file_name.split(".")[0]
    file_name = file_name.split("_")[0] #gives us the BAC ID
    add = str(str(file).split(".")[0])
    print (add)
    NewFILEPath = "%s%s" %(add,str("_reformatted.fa"))
 
    newAssembly = open(NewFILEPath, "a")
    #NewFileName = "%s%s%s"%(str(file).split(".")[0],str("_reformatted"),str(".fa"))
    for seq in SeqIO.parse(file, "fasta"):
        old_ID = str(seq.id)
        new_ID = "%s%s%s" %(str(file_name),str("_"),str(old_ID))
        new_sequence = SeqRecord(Seq(str(seq.seq)), id = new_ID, description = "")
        SeqIO.write(new_sequence, newAssembly,"fasta")
    print ("%s%s" %(str(file.split("/")[-1]), "- done!"))
    
    
    
    
