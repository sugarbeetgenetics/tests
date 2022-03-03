#!/usr/bin/python
from Bio import SeqIO
import sys
cmdargs = str(sys.argv)
for seq_record in SeqIO.parse(str(sys.argv[1]), "fasta"):
 output_line = '%s\t0\t%i\t%s\t900\t=' % \
(seq_record.id, len(seq_record), seq_record.id)
 print(output_line)
##remove\t0 from 6th line when use with FastaSummary script.
