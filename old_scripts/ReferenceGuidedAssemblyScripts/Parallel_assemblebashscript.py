import progressbar
import argparse
import os
import numpy as np
import pandas as pd
from datetime import datetime
from mpi4py import MPI
import subprocess

comm = MPI.COMM_WORLD
size = comm.Get_size() #gives of ranks in comm (number of processors to use)
rank = comm.Get_rank()

split_data = None

def split(container, count):
    #function splitting a container into equal sizes
    return [container[_i::count] for _i in range(count)]

#initial job of reading and splitting files for different process will be done in process rank=0
if rank == 0:
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-sb_f', '--file', type=str, metavar='', required=True, help='Super blocks file')
    args = parser.parse_args()
    block_file = args.file
    with open(block_file) as b_f:
        data = b_f.readlines()
        print("file read at %s" % (datetime.now()))
        print("Splitting the data!")
        split_data = split(data, size)
        print("Splitting done!")
        data_size = np.shape(data)[0]
        print(data_size)
else:
    data = None

calc = comm.scatter(split_data, root=0) #scatters the split_data to all the available processes

subprocess.call('echo "${}" | /path/to/script --args'.format(VAR))
