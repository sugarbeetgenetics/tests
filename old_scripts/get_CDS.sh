#!/bin/bash

######## install bioawk if not installed ########
#"""
#install bison to get yacc
#sudo apt install bison
#git clone https://github.com/lh3/bioawk.git
#cd bioawk/
#make
#sudo cp bioawk /usr/local/bin/
#"""
# $1 is genome.fasta
# $2 is genome.gff
# $3 is output.fasta

~/softwares/gffread/gffread/gffread -x - -g $1 $2 | \
bioawk -c fastx '{ print $name, $seq }' | \
while read line; \
do \
name=$(echo $line | cut -f 1); \
echo $line | cut -f 2 | \
awk -F "" '{ for (i = 3; i <= NF; i += 3) \
printf "%s%s", $i, (i+3>NF?"\n":FS) }' | \
awk -v name="$name" '{ print ">"name; print $1 }'; \
done \
> $3
