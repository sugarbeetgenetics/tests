#!/bin/bash
## $1 id file, $2 number of file to be made from big file, $3 number of processes
mkdir "$1_divided_file"
cd "$1_divided_file"
lines=$(< $1 wc -l)
lines_per_file=$(($lines / $2)) #$2 is number of files needed
echo "Splitting file"
#split -d -l "$lines_per_file" $1
echo "Line8 $lines"
pwd
for file in *[0-9][0-9]
do

    lines_per_small_file=$(< "$file" wc -l)

    if [[ "$lines_per_small_file" -ge $3  ]]; then
        echo Processing "$file"
        mpirun -n $3 python3 /media/pflanz/Hs1-2/Scripts/scan_samtools_depth_v2.py -mf "$file"
    else
        mpirun -n "$lines_per_small_file" python3 /media/pflanz/Hs1-2/Scripts/scan_samtools_depth_v2.py -mf "$file"
    fi
    echo Done "$file"
done

cat *Depth_Summary.txt > "$1_complete_depth_summary.txt"
# until now what this script is doing is dividing big file into small chunks based on number of line
# and then processing each file through python script to calculate depth summary for each scaffold/contig
# but what we have in Depth summary file are some scaffold data which needs to be merged
# for that merging we will use separate python script again with -tk option set to <merge> along with -f directing to
# the file needs to be processed
