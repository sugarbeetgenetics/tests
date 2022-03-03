#!/bin/bash

for i in {1..10}; do
    stage="P"$i
    echo "python3 /media/pflanz/Hs1-2/Scripts/gatk/snp_filtering.py $stage">>/media/pflanz/Hs1-2/Scripts/gatk/test.commands

done

parallel -j 12 < /media/pflanz/Hs1-2/Scripts/gatk/test.commands


