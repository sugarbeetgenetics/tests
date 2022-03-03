#####################################################
##############Author: Avneesh Kumar##################
##translate any DNA sequence in to protein sequence##
#####################################################


from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from os import path
import sys, os, glob
from progressbar import ProgressBar
pbar=ProgressBar()

#input_fasta = raw_input('Subject fasta full path:')

cmdargs = str(sys.argv)

cwd = os.getcwd()

outfile = "%s%s%s" %(str(cwd),str('/'),str('translated_seq.txt'))
with open(outfile, 'a') as w_file:
	for seq_record in SeqIO.parse(str(sys.argv[1]), "fasta"):
		template_dna = seq_record.reverse_complement()
		m_rna = (seq_record.seq).transcribe()
		protein = m_rna.translate()
		protein_seq = SeqRecord(protein)
		protein_seq.id = seq_record.id
		SeqIO.write(protein_seq, w_file,'fasta')


