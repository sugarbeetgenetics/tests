#######FINISHED

#This script is to scan samtools depth file and give average depth per scaffolds/contig while using parallel processing

import progressbar
import argparse
import os
import numpy as np
import pandas as pd
from datetime import datetime
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size() #gives of ranks in comm (number of processors to use)
rank = comm.Get_rank()

result = np.array([["#header", 0, 0]])


def split(container, count):
    #function splitting a container into equal sizes
    return [container[_i::count] for _i in range(count)]


#print("Processing Started at %s" % (datetime.now()))
file_name = str()
split_data = None

#initial job of reading and splitting files for different process will be done in process rank=0
if rank == 0:
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-mf', '--file', type=str, metavar='', required=True, help='Samtools depth file')
    args = parser.parse_args()
    depth_file = args.file
    file_name = depth_file + "_Depth_Summary.txt"
    with open(depth_file) as d_f:
        data = d_f.readlines()
        print("file read at %s" % (datetime.now()))
        print("Splitting the data!")
        split_data = split(data, size)
        print("Splitting done!")
        data_size = np.shape(data)[0]
        print(data_size)
else:
    data = None

calc = comm.scatter(split_data, root=0) #scatters the split_data to all the available processes


for line in progressbar.progressbar(calc):
    line_split = np.empty((1, 3))
    if line == "\n":
        pass
    else:
        line_split = np.array([str(line).split(sep="\t")])
        line_split[0, 2] = str(line_split[0, 2]).split(sep="\n")[0]
    #check for scaffold, if it is present update it, if not then add new row
    if line_split[:, 0] in result[:, 0]:
        scaf_data = result[np.where(result[:, 0] == line_split[0, 0])]
        scaf_depth = int(scaf_data[0, 1]) + int(line_split[0, 2])
        scaf_freq = int(scaf_data[0, 2]) + 1
        result[np.where(result[:, 0] == line_split[0, 0])[0].tolist(), 0] = line_split[0, 0]
        result[np.where(result[:, 0] == line_split[0, 0])[0].tolist(), 1] = scaf_depth
        result[np.where(result[:, 0] == line_split[0, 0])[0].tolist(), 2] = scaf_freq
    else:
        result = np.vstack((result, line_split))

results = np.array(comm.gather(result, root=0))

final_result = np.array([["#header", 0, 0]])

if rank == 0:
    print("final calculation at %s" % (datetime.now()))

    reduced_dim_results = results.reshape((np.shape(results)[0]*np.shape(results)[1], np.shape(results)[2]))
    for row in reduced_dim_results:
        if row[0] in final_result[:, 0]:
            final_result[np.where(final_result[:, 0] == row[0])[0].tolist(), 0] = row[0]
            final_result[np.where(final_result[:, 0] == row[0])[0].tolist(), 1] = int(final_result[np.where(
                final_result[:, 0] == row[0])[0].tolist(), 1]) + int(row[1])
            final_result[np.where(final_result[:, 0] == row[0])[0].tolist(), 2] = int(final_result[np.where(
                final_result[:, 0] == row[0])[0].tolist(), 2]) + int(row[2])
        else:
            final_result = np.vstack((final_result, row))
#    print((final_result))

    cwd = os.getcwd()
    df = pd.DataFrame(final_result)
    with open(str(file_name), 'w') as out_file:
        df.to_csv(out_file, sep="\t", header=False, index=False)

#    print(datetime.now())
