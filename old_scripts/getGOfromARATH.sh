#!/bin/bash

awk '{print $1}' $1 | grep -w -f - /media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/annotation/final_annotated_v4.3.5_sorted.gff3 | awk -F"\t" '{split($9,a,"|"); if($3=="mRNA") print a[4] }' | grep -w -f - /media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/annotation/ARATH_3702_idmapping.dat.Gene_OrderedLocusName | awk '{split($3,a,"_"); if(a[2]=="") print $3; else print a[2]}' > $1_AtIDs.txt
echo "Done $1"
