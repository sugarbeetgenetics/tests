#!/bin/bash

## $1 windows file, $2 number of file to be made from big file, $3 number of processes, $4 depth file
mkdir "$PWD/divided_window_file"
cd "$PWD/divided_window_file"
lines=$(< $1 wc -l)
lines_per_file=$(($lines / $2)) #$2 is number of files needed
echo "Splitting file"
split -d -n l/"$lines_per_file" $1
echo "$lines"
pwd
for file in *[0-9][0-9]
do

    lines_per_small_file=$(< "$file" wc -l)

    if [[ "$lines_per_small_file" -ge $3  ]]; then
        echo Processing "$file"
        mpirun -n $3 python3 /media/pflanz/Hs1-2/Scripts/Find_Translocation/windowed_Coverage.py -wf "$file" -df $4
    else
        mpirun -n "$lines_per_small_file" python3 /media/pflanz/Hs1-2/Scripts/Find_Translocation/windowed_Coverage.py -wf "$file" -df $4
    fi
    echo Done "$file"
done

cat *_Windowed_Depth_Summary.txt > "$PWD/complete_window_depth.txt"