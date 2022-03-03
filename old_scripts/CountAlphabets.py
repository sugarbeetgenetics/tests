#This script count #N in a sequence/s

from Bio import SeqIO
import argparse
import os
import csv

parser = argparse.ArgumentParser(description='Count #N in given sequence and the length of Seq')
parser.add_argument('-s', '--sequence', type=str, metavar='', required=True, help='single fasta or multifasta')
args = parser.parse_args()

seq_file = args.sequence
lengths = []
counts = []
ids = []
percentNs = []

for seq in SeqIO.parse(seq_file, "fasta"):
    i = 0
    count = sum(str(seq.seq).count(x) for x in ['N', 'n'])
    if count != 0:
        counts.append(count)
        length = len(seq.seq)
        lengths.append(length)
        ids.append(seq.id)
        percentN = ((float(count)/float(length))*100)
        percentNs.append(int(percentN))

cwd = os.getcwd()

data = zip(ids, lengths, counts, percentNs)

with open(str(cwd+'/'+'counts.txt'), 'w') as out_file:
    for x in data:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerows(data)
