#!/bin/bash

CWD=$(pwd)
for file in "$CWD"/*[12]
do
  fname=${file##*/} #This gives your base filename.
  fpath=${file%/*} # Your dir
  dname=${fpath##*/}
  echo $fpath/$fname
  cd $fpath/$fname
  ls
  cp test1.txt ${fpath}/${fname}/${fname}.test1.txt
  cd ../
done