import pandas
import csv

#############################################
########### table processing ################
# remove heading
# remove tab at line ends (replace \t\n with \n)
# write "ID" for first column
#############################################

infile = "/media/pflanz/BackUps/Brassicanapus/DataFromChina/HaplotypeAnalysis/genoploteR/Brapa_napus.raw.filter.snp.anno.plinl.chrA10_13349572-13360649_downstream_targets_of_FLC.vcf.flpjk_Default View_565x36.dat"
outfile = "/media/pflanz/Avneesh/Haplotypes_output.csv"

variantFile = pandas.read_csv(infile, sep="\t", comment = "#")

########################################################
########################################################
# function to filter for missing values
def filter_missing(all_snps, percent_missing_allowed):
    missing = 0
    for item in all_snps:
        if item == "-":
            missing += 1
    missing_percent = (missing*100)/len(all_snps)
    #return missing_percent
    if missing_percent <= percent_missing_allowed:
        return False
    elif missing_percent > percent_missing_allowed:
        return True
#######################################################
haplotypes = {}
haplotypes_count = 0
ids_included = []
filtered_missing = 0
percent_missing_allowed = 5
max_missmatch_allowed = 1

i=0
while i < (variantFile.shape[0]): #(variantFile.shape[0])
    #print("i:" + str(i))
    all_snps = list(variantFile.iloc[i,1:])
    #print("all_snps:"  + str(all_snps))
    j = 0
    #print("line1")
    # filter for missing
    if filter_missing(all_snps,percent_missing_allowed):
        filtered_missing += 1
        i += 1 #it is needed here because from this point loop will go to the top, so there is need to update 
        # i at this stage
        #print("line2")
        continue #skip further parts of this loop if this row needs to be filtered

    if variantFile.iloc[i,0] not in ids_included: # to run only for ids which are not in any haplotype
        #print("line3")
        same_as = [variantFile.iloc[i,0]]
        #ids_included.append(variantFile.iloc[i,0])

        while j <  variantFile.shape[0]: #variantFile.shape[0]
            #print("j:" + str(j))
            #print("line4")
#            if filter_missing(all_snps,10):
#                filtered_missing += 1
                #print("line2")
#                continue #skip further parts of this loop if this row needs to be filtered

            if i!=j:
                all_snps2 = list(variantFile.iloc[j,1:])
                #print("all_snps2:" + str(all_snps2))
                #print("line5")
                if filter_missing(all_snps2, percent_missing_allowed):
                    #print("line6")
                    j += 1
                    continue #skip further parts of this loop if this row needs to be filtered
                if variantFile.iloc[j,0] not in ids_included: # to run only for ids which are not in any haplotype
                    #print("line3")
                    #same_as = [variantFile.iloc[i,0]]
                    #ids_included.append(variantFile.iloc[i,0])

                    match_count = 0
                    missmatch_count = 0
                    k = 0
                    while k < len(all_snps):
                        #print("line7")
                        if all_snps[k]==all_snps2[k] and (all_snps[k]=="-" or all_snps2[k]=="-"):
                            match_count += 1
                        elif (all_snps[k]!=all_snps2[k] and all_snps2[k]!="-"):
                            missmatch_count += 1
                        k += 1
                    if missmatch_count < max_missmatch_allowed : #missmatch_count*100/len(all_snps)
                        #print("line8")
                        same_as.append(variantFile.iloc[j,0])
               
            j += 1
        ids_included.extend(same_as[:])
        haplotypes[haplotypes_count] = same_as
        haplotypes_count += 1
            
    i += 1
    
# write down dictionary (haplotypes) to a file
dict = haplotypes

w = csv.writer(open(outfile, "w"))
for key, val in dict.items():
    w.writerow([key, val])


