
##### This script find gaps in sequence

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, os, glob
import csv

#Assemblies = glob.glob("/media/avneesh/AneeshHDDfat/AssembledScaffolds/*")

args = sys.argv

file = args[1]
product_size = int(args[2])

all_gaps_start = []
all_gaps_end = []

gap_file = str(file) + str("_gaps_info.csv")

for seq in SeqIO.parse(file, "fasta"):
    #print("seq length:%s"%len(seq.seq))
    base = 0
    seq_end = "no"
    new_sub_number = 0
    for a in range(base, len(seq.seq),1):
        if seq.seq[a] == "N":
            if (a>0 and seq.seq[a-1] == "N"):
                #print("line 32")
                if (a>0 and seq.seq[a-1] == "N" and a==len(seq.seq)-1):
                    gap_end = a
                    #print("line 45:%s"%gap_end)
                    all_gaps_end.append(gap_end)
            elif (a>0 and seq.seq[a-1] != "N"):
                gap_start = a
                #print("line 34:%s"%gap_start)
                all_gaps_start.append(gap_start)
        else:
            if (a>0 and seq.seq[a-1] == "N"):
                gap_end = a-1
                #print("line 41:%s"%gap_end)
                all_gaps_end.append(gap_end)                
            elif (a>0 and seq.seq[a-1] != "N"):
                #print("line 48")
                pass
        #print("a:%s"%a)

#print("new_all_start: %s"%len(all_gaps_start))
#print("new_all_end: %s"%len(all_gaps_end))

with open(gap_file, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(zip(all_gaps_start, all_gaps_end))


####regions to be included for primer designing

#first we expand the gap by 100bp to both the sides
new_all_start = []
new_all_end = []

for start in all_gaps_start:
    if start>100:
        new_start = start - 100
    else:
        new_start = 0
    new_all_start.append(new_start)    
for seq in SeqIO.parse(file, "fasta"):
    seq_length=len(seq.seq)
    for end in all_gaps_end:
        if end+100<len(seq.seq):
            new_end = end+100
        else:
            new_end = len(seq.seq)
        new_all_end.append(new_end)

stretched_gaps = str(file) + str("_stretchedgaps_info.csv")
with open(stretched_gaps, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(zip(new_all_start, new_all_end))

#print("new_all_start: %s"%len(new_all_start))
#print("new_all_end: %s"%len(new_all_end))

#define sequence_included_region which should be max 2000bp including gap size
all_amplification_start = []
all_amplification_end = []

for i in range(0, len(new_all_start),1):
    start = int(new_all_start[i])
    end = int(new_all_end[i])
    if (end-start)<product_size:
        to_be_amplified = product_size-(end-start)
    expansion = int(round(to_be_amplified/2))
    #print(expansion)
    if (start-expansion > 0):
        amplification_start = start-expansion
    else:
        amplification_start = 0
    if (end+expansion < seq_length):
        amplification_end = end+expansion
    else:
        amplification_end = seq_length

    all_amplification_start.append(amplification_start)
    all_amplification_end.append(amplification_end)
    
amplificationregion = str(file) + str("_region_info.csv")

with open(amplificationregion, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(zip(all_amplification_start, all_amplification_end))


