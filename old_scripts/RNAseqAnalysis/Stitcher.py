#!/usr/bin/python

#Stitcher.py
#can stitch .gff file with FPKM value file

__author__ = 'Avneesh Kumar'

import argparse
from itertools import chain
import csv, sys

parser = argparse.ArgumentParser(description='This script can stitch .gff file with FPKM value file. --Avneesh Kumar')
parser.add_argument('-gff','--gff_file',help='Input gff file path', required=True)
parser.add_argument('-exp','--FPKM_file',help='Input FPKM file path', required=True)
parser.add_argument('-fra','--fraction',help='Fraction of gene length to increase gene size (to incorporate UTRs) (default=20)', required=False, default=20, type=int)
parser.add_argument('-min','--minimum_FPKM',help='minimum FPKM default=0', required=False, default=0, type=int)

args = parser.parse_args()

#open gff
#gff = open(args.gff_file, 'r')

tralo = ['scaffold4317','scaffold1723','scaffold11953','C5063194','scaffold33615','scaffold12163','scaffold7599','scaffold19719','scaffold15066','scaffold25850','scaffold39869','scaffold19738','scaffold26666','scaffold5217','scaffold13756','scaffold8090','scaffold32842']
f = args.fraction
min = args.minimum_FPKM
#read files in loop
i=0
with open(args.gff_file) as gff:
    for gff_line in gff:
        i=i+1
        tr = 0
        scaffold = gff_line.split("\t")[0]
        start = int(gff_line.split("\t")[3])
        end = int(gff_line.split("\t")[4])
        length = abs(end-start)
        strand = gff_line.split("\t")[6]
        name = gff_line.split(";")[-1]
        name = name[0:-1] #to remove new line character which is coming along from input filer
        #expand the gene region
        if length>500:
            modification = int(length/f)
        else:
            modification = int(2*length)

        start_out = start - modification
        start_in = start + modification
        end_out = end + modification
        end_in = end - modification

        if scaffold in tralo:
            tr = 1
        else:
            tr = 0
        with open(args.FPKM_file) as fpkm:
            for (line) in fpkm:
                if scaffold in line:
                    position = line.split("\t")[6]
                    coordinates = position.split(":")[1]
                    fpkm_start = int(coordinates.split("-")[0])
                    fpkm_end = int(coordinates.split("-")[1])
                    fpkm_value = float(line.split("\t")[9])
                    #print ("%s\t%d\t%d\t%s\t%d\t%d" %(scaffold,start,end,position,fpkm_start,fpkm_end))
                    if (start_out < fpkm_start < end_out ) and fpkm_value > min : #or end_out > fpkm_end > end_in
                            sys.stdout.write("%s\t%s\t%d\t%d\t%f\t%s\t%d\n" %(scaffold,coordinates,start,end,fpkm_value,name,tr))
