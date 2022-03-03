#!/bin/bash

# $1 and $2 are fastq read file and $3 is read IDs list.
fastq_out1=/media/pflanz/Avneesh/Test_Bash/R1.fastq
fastq_out2=/media/pflanz/Avneesh/Test_Bash/R2.fastq
out=/media/pflanz/Avneesh/Test_Bash

/home/pflanz/seqtk/seqtk subseq $1 $3 1>$fastq_out1
/home/pflanz/seqtk/seqtk subseq $2 $3 1>$fastq_out2
~/SPAdes-3.12.0-Linux/bin/spades.py -t 5 -m 100 -k 39 --pe1-1 $fastq_out1 --pe1-2 $fastq_out2 -o $out
