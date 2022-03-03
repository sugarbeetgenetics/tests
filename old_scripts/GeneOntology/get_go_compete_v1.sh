#!/bin/bash
cwd=$(pwd)
awk '{if($3=="mRNA") print $0}' $1 | awk '{n=split($9,a,",|;"); for (i=1; i < n; ++i) if(a[i] ~ "^GO") print $1":"$4"-"$5"\t" a[i]}' > $(pwd)/outfile_go.txt

python3.6 /media/pflanz/Hs1-2/Scripts/GeneOntology/get_annotation.py -l $(pwd)/outfile_go.txt

python3.6 /media/pflanz/Hs1-2/Scripts/AdvancedPlots/wordcloud_generate.py -d  $(pwd)/v1_go_Annotation.out
