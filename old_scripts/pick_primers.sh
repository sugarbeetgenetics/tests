#!/bin/sh
read gene_list
read primer3_config_file
while read line; do    
    awk -F "\t" '{}' $line    
done < file.txt
