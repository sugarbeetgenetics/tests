#!/bin/bash

# $1 id superblockfile and $2 id bam file, $3 is fatsqR1, $4 is fastqR2.

#out=
sp_block_File=$1
#mkdir "$1_Extract_Read"
#cd "$1_Extract_Read"
while read LINE; do
    while IFS=$'\t' read -r -a ADDR; do
        #echo "${ADDR[0]}:${ADDR[1]}-${ADDR[2]}"
        samtools view $2 "${ADDR[0]}:${ADDR[1]}-${ADDR[2]}" |
        {
        while IFS=$'\t' read -r -a SAMDATA; do
           touch "${ADDR[0]}:${ADDR[1]}-${ADDR[2]}_IDS.lst"
           echo "${SAMDATA}" >>"${ADDR[0]}:${ADDR[1]}-${ADDR[2]}_IDS.lst"
        done
        }
        #echo $fastq_out1 $fastq_out2 $READDATA
        fastq_out1="${ADDR[0]}:${ADDR[1]}-${ADDR[2]}_R1.fastq"
        fastq_out2="${ADDR[0]}:${ADDR[1]}-${ADDR[2]}_R2.fastq"
        /home/pflanz/seqtk/seqtk subseq $3 "${ADDR[0]}:${ADDR[1]}-${ADDR[2]}_IDS.lst" 1>$fastq_out1
        /home/pflanz/seqtk/seqtk subseq $4 "${ADDR[0]}:${ADDR[1]}-${ADDR[2]}_IDS.lst" 1>$fastq_out2
        #/media/pflanz/Hs1-2/Scripts/ReferenceGuidedAssemblyScripts/fastq-pair/build/fastq_pair $fastq_out1 $fastq_out2
        ~/SPAdes-3.12.0-Linux/bin/spades.py -t 5 -m 100 -k 39 --pe1-1 $fastq_out1 --pe1-2 $fastq_out2 -o "./${ADDR[0]}:${ADDR[1]}-${ADDR[2]}_out" >/dev/null 2>&1
    done <<< "$LINE"
done < $sp_block_File


