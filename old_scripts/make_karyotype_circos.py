#!/usr/bin/python

import random
from Bio import SeqIO
import sys
cmdargs = str(sys.argv)
###use next two lines for random colors
#color_file = open("/media/planz/Avneesh/Beet_translocation/ProperDataSet/CIRCOS/colors.txt","r") 
#color_ids = color_file.readlines()

color_ids2 = ['black','red'] 

for seq_record in SeqIO.parse(str(sys.argv[1]), "fasta"):
	#color_id = random.choice(color_ids)
	output_line = 'chr\t-\t%s\t%s\t0\t%i\tred' % \
	(seq_record.id,seq_record.id, len(seq_record))  # use ''', color_id.split('\n')[0]''' to use random colors
	print(output_line)
##add color manually in the output file.
#replace red with %s to use random colors
