#this scripts make sense out of primer blast data

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys, os
import csv, re

args = sys.argv
primer_fasta = args[1]
blastn_file = args[2]
fasta_id = args[3]
primers_info = args[4]
#print (primers_info)

#blastn_file = os.environ["primers"]
#primer_fasta = int(os.environ["blastn"])

### find NGG sites ###



#############get seq ids of bad primers############

bad_primer_id = []
print ("Started")
for seq in SeqIO.parse(primer_fasta, "fasta"):
    #print (seq)
    with open(blastn_file, "r") as search:
        for line in search:
            line = line.rstrip()
            #print (seq.id)
            if line.startswith(seq.id): 
                #print (line) 
                #print (line.split("\t")[2])
                if int(float(line.split("\t")[2])) == 100:
                    #print (line.split("\t")[2])
                    if int(float(line.split("\t")[3])) == len(seq.seq):
                        #print (line)
                        #print ("Good")
                        if line.split("\t")[1] != fasta_id:
                            bad_primer_id.append(int(float(line.split("\t")[0])))

#print (bad_primer_id)


######################
cwd=os.getcwd()

###### keep only good primers #################

for seq in SeqIO.parse(primer_fasta, "fasta"):
    #print (seq.id)
    if int(seq.id) not in bad_primer_id:
        with open("useful_primers.fasta","a") as useful_primers:
            new_ID = seq.id
            new_sequence = SeqRecord(Seq(str(seq.seq)), id = new_ID, description = "")
            SeqIO.write(new_sequence, useful_primers,"fasta")

i=0
for seq in SeqIO.parse(primer_fasta, "fasta"):
    #print (seq.id)
    if int(seq.id) not in bad_primer_id:
        #path = str(cwd+seq.id+"_PRIMER_INFO")
        print (primers_info)
        #f=open((cwd+"/"+primers_info+"_PRIMER_INFO"), 'r')    
        #primers_info=f.read()
        for line in primers_info:
            print (line)
            sequence=str(seq.seq)
            print (sequence)
            #if re.search(re.escape(sequence), line, re.I):
            #    with open("useful_primers_info.fasta","a") as useful_primers_info:
            #        useful_primers_info.write("%s\n"%line)
            #i +=1
#(r"^"+re.escape(seq.seq)+r"?")
