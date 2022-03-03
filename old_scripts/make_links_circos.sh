#!/bin/bash

echo "Enter the file fath(BLAST format7)"
read -p 'file: ' blast_file
read -p 'minimum %match: ' match
read -p 'minimum length of match: ' length
echo "Creating links data. Have a cup of coffee... !"

sed '/#/d' $blast_file | awk '($3>"$match") && ($4>"$length") {print $1,"\t",$7,"\t",$8,"\t",$2,"\t",$9,"\t",$10}' > ./links.txt
