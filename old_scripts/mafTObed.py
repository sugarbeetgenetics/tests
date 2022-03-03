#!/usr/bin/python

'''
Author: Avneesh Kumar
At Plant Breeding Institute, CAU Kiel

make simple bed file from simple maf(like) file.")
'''

from sys import argv
import csv



with open(new_file, 'wb') as outf:
	writer = csv.writer(outf,delimiter='\t')
	writer.writerows(izip(ids,win_start,win_end))
			
		
print('Finished!')
