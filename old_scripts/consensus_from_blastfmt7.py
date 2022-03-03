###############################################################
#### Take blast output (fmt 7) file and extend the         ####
#### sequence for overlapping scaffolds.                   ####
###############################################################

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from progressbar import ProgressBar
pbar=ProgressBar()

#blast_output = raw_input('Blast output file:')
#fasta_file = raw_input('Fasta fiel:')

sequence_fasta_file = '/media/planz/Hs1-2/Beet_translocation/ProperDataSet/BLASTdb/P_procumbens-1.0.soap.1.fa'
query_fasta_file = '/media/planz/Hs1-2/Beet_translocation/ProperDataSet/BLAST/A_minus_B/ScaffoldSeq_OnlyTR520.fasta'
blast_output = '/media/planz/Hs1-2/Beet_translocation/ProperDataSet/BLAST/A_minus_B/ScaffoldSeq_OnlyTR520_vs_procumbens_cleaned_awk.txt'
n = 10 #leave n nucleotides from the edges

bfile = open(blast_output)
bfile_lines = bfile.readlines()

q_id = []
s_id = []
q_start = []
q_end = []
s_start = []
s_end = []

for yz in bfile_lines:
    q_id.append(yz.strip().split('\t')[0])
    s_id.append(yz.strip().split('\t')[1])
    q_start.append(yz.strip().split('\t')[2])
    q_end.append(yz.strip().split('\t')[3])
    s_start.append(yz.strip().split('\t')[4])
    s_end.append(yz.strip().split('\t')[5])
    
print ('Blast output read!')
q_record = SeqIO.to_dict(SeqIO.parse(query_fasta_file, "fasta"))
s_record = SeqIO.to_dict(SeqIO.parse(sequence_fasta_file, "fasta"))

print ('Fasta read!')

def useful_scaffold(q_id,s_id,q_start,q_end,s_start,s_end):
    global q_record, s_record, n
    status = 0
    q_len = len(q_record[q_id].seq)
    s_len = len(s_record[s_id].seq)
    if s_start<s_end:
        if (q_id==s_id and q_start==s_start and q_end==s_end):
            status = 0
        elif (q_start>n and q_end==q_len and s_start<=n and s_end==s_len):
            status = 0
        elif (q_start>n and q_end==q_len and s_start<=n and s_end<s_len):
            status = 1
        elif (q_start<=n and q_end<q_len and s_start<=n and s_end==s_len):
            status = 0
        elif (q_start<=n and q_end<q_len and s_start>1 and s_end==s_len):
            status = 1
        elif (q_start>1 and q_end<q_len and s_start<=n and s_end==s_len):
            status = 0
            
    elif s_start>s_end:
        if (q_start<=n and s_start==s_len and q_end==q_len and s_end<=n):
            status = 0
        elif (q_start<=n and s_start==s_len and q_end<q_len and s_end<=n):
            status = 0
        elif (q_start<=n and s_start<s_len and q_end<q_len and s_end<=n):
            status = 1
        elif (q_start>1 and s_start==s_len and q_end==q_len and s_end<=n):
            status = 0        
        elif (q_start>1 and s_start==s_len and q_end==q_len and s_end>1):
            status = 1
        elif (q_start>1 and s_start==s_len and q_end<q_len and s_end<=n):
            status = 0
    return(status)

 
scaffold_status = []   
for i in pbar(range(0,len(q_id))):
    if (useful_scaffold(q_id[i],s_id[i],q_start[i],q_end[i],s_start[i],s_end[i])==1):
        print(scaffold_status)


with open('/media/planz/Hs1-2/Beet_translocation/ProperDataSet/BLAST/A_minus_B/scaffold_status.txt',"w") as fileout:
    if (useful_scaffold(q_id[i],s_id[i],q_start[i],q_end[i],s_start[i],s_end[i])==1):    
        for q_id[i],s_id[i],q_start[i],q_end[i],s_start[i],s_end[i] in zip(q_id[i],s_id[i],q_start[i],q_end[i],s_start[i],s_end[i]):
            fileout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(q_id[i],s_id[i],q_start[i],q_end[i],s_start[i],s_end[i]))

'''
def make_new_seq(q_id,s_id,q_start,q_end,s_start,s_end):
    global record
    q_len = len(record[q_id].seq)
    s_len = len(record[s_id].seq)
    new_seq = ''
    if s_start < s_end:
        new_seq = new_seq + record[q_id] + (record[s_id].seq)[s_end+1:-1]
    elif s_start > s_end:
        new_seq = new_seq + (record[s_id].seq)[1:s_start] + record[q_id].seq
    return(str(new_seq))

cwd = os.getcwd()
output_file = "%s%s%s" %(str(cwd),str("/new_combined_scaffolds"),str(".fa"))

 
with open(output_file, 'a') as w_file:
    for i in pbar(range(0,len(q_id))):
        if (useful_scaffold(q_id[i],s_id[i],q_start[i],q_end[i],s_start[i],s_end[i])==1):
            seq = Seq(make_new_seq(q_id[i],s_id[i],q_start[i],q_end[i],s_start[i],s_end[i]))    
            SeqIO.write(SeqRecord(seq, id=(q_id[i]+'_'+s_id[i])), w_file,'fasta')


'''


