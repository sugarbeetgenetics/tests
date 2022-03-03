'''
Author: Avneesh Kumar
At Plant Breeding Institute, CAU Kiel

make simple bed file for multi-fasta file. make bins for specified size. change window size below in line 22")
'''

from os import path
import sys, os, glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
from itertools import izip

from progressbar import ProgressBar
pbar=ProgressBar()

args = sys.argv

infile = '/media/planz/Hs1-2/Beet_translocation/ProperDataSet/BLASTdb/RefBeet-1.2.fna' ####insert file path
win_size = 100 ####insert windows size


#make dictionary of scaffolds and their size

seq_lengths = {}
for record in SeqIO.parse(infile, "fasta"):
	seq_lengths[record.id] = len(record.seq)
	
print('Read fasta!')	
	

win_start = []
win_end = []
ids = []

for seq_id,seq_len in seq_lengths.iteritems():
	print(seq_id)
	if seq_len%2 == 0:
		windows = seq_len/win_size
	else:
		windows = (seq_len/win_size)+1
	
	st = 0
	end = win_size
	win = 0
	for win in range(0,windows):
		ids.append((seq_id+' length='+str(seq_len)))
		win_start.append(st)
		win_end.append(end)
		
		if (end+win_size) < seq_len:
			st = end+1
			end = end + win_size
		else:
			print('Exiting')
			break

print ('Success!')
		
new_file = "%s%s%s%s%s" %(str(infile).split('.')[0],str('_'),str('bins'),str('.'),str('bed'))

with open(new_file, 'wb') as outf:
	writer = csv.writer(outf,delimiter='\t')
	writer.writerows(izip(ids,win_start,win_end))
			
		
print('Finished!')
