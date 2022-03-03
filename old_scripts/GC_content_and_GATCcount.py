###by Avneesh Kumar###
######usage####
#after first argument on command line should be the fasta file, second should be "plot" (if you want scatterplot),
#and third should be to specify if you need sequences with length less than 1000 or more than 1000 (and in case you
# have another length to distinguish sequences into two groups, just change "1000" in line 17 and 25 with your number)

from Bio import SeqUtils
from Bio import SeqIO
import re
import sys
import csv
import os
import progressbar
cmdargs = str(sys.argv)


id = []
length = []
gc = []
gatc = []
to_find = "GATC"
pattern = re.compile(to_find, re.IGNORECASE)

for record in progressbar.progressbar(SeqIO.parse(str(sys.argv[1]),"fasta")):
    i = 0
    id.append(record.id)
    seq = str(record.seq)
    length.append(len(str(record.seq)))
    gc.append(SeqUtils.GC(record.seq))
    for match in pattern.finditer(seq):
        i = i+1
    gatc.append(i)

cwd = os.getcwd()
with open(str(cwd+'/'+'GATC_GC_counts.txt'), 'w') as out_file:
    writer = csv.writer(out_file, delimiter='\t')
    writer.writerows(zip(id,length,gc,gatc))