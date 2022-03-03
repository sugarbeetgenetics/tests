###by Avneesh Kumar###
######usage####
#after first argument on command line should be the fasta file, second should be "plot" (if you want scatterplot),
#and third should be to specify if you need sequences with length less than 1000 or more than 1000 (and in case you
# have another length to distinguish sequences into two groups, just change "1000" in line 17 and 25 with your number)

from Bio import SeqUtils
from Bio import SeqIO
import matplotlib.pyplot as plt
import sys
cmdargs = str(sys.argv)


id = []
gc = []
gatc = []

for record in SeqIO.parse(str(sys.argv[1]),"fasta"):
    if(str(sys.argv[3])=="less" and len(record.seq)<1000):
        output_line = '%s\t%i' % \
        (record.id, SeqUtils.GC(record.seq))
        id.append(record.id)
        gc.append(SeqUtils.GC(record.seq))
        print (output_line)


    elif(str(sys.argv[3])=="more" and len(record.seq)>1000):
        output_line = '%s\t%i' % \
        (record.id, SeqUtils.GC(record.seq))
        id.append(record.id)
        gc.append(SeqUtils.GC(record.seq))
        print (output_line)


if(str(sys.argv[2])=="plot"):
    plt.scatter(gc,id)
    plt.show()