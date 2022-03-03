file = "/media/pflanz/Avneesh/Brassicanapus/gene_exp_significant_geneID_mRNA_PFam_finallist.txt"
with open(file) as testFile:
	newlist =[]
	for line in testFile:
		newline = []
		#print (line.rstrip())
		for item in line[0:-2].rstrip().split("-"):
			if (item in newline):
				print("")
			else:
				newline.append(item)
		newlist.append(newline)
        
with open('./uniqueValues.txt', 'w') as f:
    for item in newlist:
        f.write("%s\n" % item)
