import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
sequence_fasta_file = '/media/planz/Hs1-2/Beet_translocation/ProperDataSet/BLAST/A_minus_B/ScaffoldsTR520_notinRefBeetv1.2_inProcumbens.fasta'
sizes = [len(rec) for rec in SeqIO.parse(sequence_fasta_file,"fasta")]

#sequences in their decreasing order of size

bp_10k = []
for rec in SeqIO.parse(sequence_fasta_file,"fasta"):
    if len(rec.seq)>10000:
        bp_10k.append(rec.id)

import pylab
pylab.hist(sizes, bins=300)
pylab.title("%i Sequences\nlengths %i to %i" \
    %(len(sizes),min(sizes),max(sizes)))
pylab.xlabel("Sequence length (bp)")
pylab.ylabel("Count")
pylab.xlim([0,(max(sizes)+1000)])
#pylab.show()

print(len(bp_10k))

with open('/media/planz/Hs1-2/Beet_translocation/ProperDataSet/BLAST/A_minus_B/ScaffoldsTR520_notinRefBeetv1.2_inProcumbens_grtr_10kbp.txt','w') as file_out:
    for i in bp_10k:
        file_out.write("{}\n".format(i))