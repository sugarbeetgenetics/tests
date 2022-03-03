###################################################################
#### If A and B are two datasets for sequence scaffolds,       ####
#### this script will give scaffolds present in A but not in B.####
###################################################################


import shlex, os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing.dummy import Pool as ThreadPool

#DataA = raw_input('Dataset A:')
#DataB = raw_input('Dataset B:')
#DataA is Database
#DataB is Query

DataA = '/media/planz/Hs1-2/Beet_translocation/ProperDataSet/BLAST/SeqScaffolds/AllScaffolds.fa'
DataB = '/media/planz/Hs1-2/Beet_translocation/ProperDataSet/BLASTdb/P_procumbens-1.0.soap.1.fa'


#define funciton for blast
def run_blast():
    print ('Searching for data A scaffolds in data B.... \n')
    
    global DataB
    global DataA
    list_scaffold_onlyA =  []
    blastn = shlex.split('blastn -word_size 100 -perc_identity 98 -outfmt 7 -num_threads 8 -max_target_seqs 1 -db')
    blastn_query = blastn + [DataB]

    sed_clean = ['sed', '/^#/ d']#to clean hash tag lines
    awk_clean = ['awk','{print $1}']
    cwd = os.getcwd()
    blast_out = "%s%s%s%s" %(str(cwd),"/",str('blast_out'),str(".txt"))
    clean_blast_out = "%s%s%s%s" %(str(cwd),"/",str('clean_blast_out'),str(".txt"))

    blast_command = blastn_query + ['-query'] + [DataA] + ['-out'] + [blast_out]
    blast_output = subprocess.call([str(x) for x in blast_command])
    SED_final = SEDclean + [blast_out]
    SED_output = subprocess.check_output([str(y) for y in SED_final])
    awk = awk_clean + [SED_output]
    
    if SED_output != '':
        list_scaffold_onlyA.append('Yes\n')
        print('Yes\n')
        os.remove(blast_out)
        os.remove(temp_seq)
    else:
        list_scaffold_onlyA.append('No\n')
        print('No\n')
        os.remove(blast_out)
        os.remove(temp_seq)
        break
    return(blast_output)

cwd = os.getcwd()
Scaffold_OnlyA = "%s%s%s%s" %(str(cwd),"/",str('Scaffold_OnlyA'),str(".txt"))

with open(Scaffold_OnlyA,"w") as onlyA:
    onlyA.writelines(run_blast())

