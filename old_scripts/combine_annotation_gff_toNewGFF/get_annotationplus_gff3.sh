#!/bin/bash


# 1. tmhmm.out must be filtered using following shell command (last column tells where the domain is)
awk '{split($6,a,"="); if(a[2]~/^[io][0-9]/) print $0}' tmhmm.out > tmhmm_filtered.out
# 2. PFAM out file filter using following command ($14 is score)
awk '{if($14>20) print $1"\t"$2"\t"$4}' ../trinotate/data_from_server/protein_CDS_PFAM_filtered.out > ../trinotate/data_from_server/protein_CDS_PFAM_filtered_1_2_4.out
# 3. gff modification in order to `get gene ids in the last column`
grep -v '#' augustus.hints.gff | awk -F"\t" '{if($3=="transcript") print $0 "\t" $9 ; else print $0 "\t" "-"}' > gff_modified_temp.gff
#####grep -v '#' augustus.hints.gff | awk -F"\t" '{split($9,a,";"); split(a[1],b,"\""); if(($3=="transcript") || ($3=="gene")) print $0 "\t" "-" ; else print $0 "\t" b[2]}' > gtf_modified_temp.gff
# 4. signalp format
awk '{print $1 "\t" "SignalP"}' signalp.out > signalp_filtered.out
# 5. blastp 
awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' ../../trinotate/data_from_server/blastp.outfmt6 > ../../trinotate/data_from_server/blastp_filtered.outfmt6

#####new blast , saved in a folder named new_blast in side annotation folder
#make blast db
~/ncbi-blast-2.7.1+/bin/makeblastdb -in ../../../../databases/all_combined.fa -dbtype prot
#blast our prot seq against above generated database
~/ncbi-blast-2.7.1+/bin/blastp -query ../../trinotate/protein_eachCDS.fasta -db ../../../../databases/all_combined.fa -num_threads 12 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > proteinCDS_At_Bna_Bol_Brap_Sina.outfmt6
##use this blast result in python script to combine with gff, instead of blast against whole uniprot database.

python3 /media/pflanz/Hs1-2/Scripts/combine_annotation_gff_toNewGFF/combine_all_toGff.py -bp /media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/trinotate/data_from_server/blastp_filtered.outfmt6 -pf /media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/trinotate/data_from_server/protein_CDS_PFAM_filtered_1_2_4.out -tm /media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/trinotate/data_from_server/tmhmm_filtered.out -sp /media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/trinotate/data_from_server/signalp_filtered.out -gff /media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/braker/Bna_Express/gff_modified_temp.gff
#for new blast (only BrNapus, BrRapa, BrOla, Sinapis, Ath)
python3 /media/pflanz/Hs1-2/Scripts/combine_annotation_gff_toNewGFF/combine_all_toGff.py -bp /media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/annotation/new_blast/proteinCDS_At_Bna_Bol_Brap_Sina_filtered.outfmt6 -pf /media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/trinotate/data_from_server/protein_CDS_PFAM_filtered_1_2_4.out -tm /media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/trinotate/data_from_server/tmhmm_filtered.out -sp /media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/trinotate/data_from_server/signalp_filtered.out -gff /media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/braker/Bna_Express/gff_modified_temp.gff

#### to get final annotated dataset
#not used awk -F"\t" '{if($3=="transcript") print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "transcript_id" " \"" $9 "\"; pfam_id \"" $11 "; UniProt_id \"" $14 "^" $15 "^" $16 "\" ;tmhmm \"" $18 "\""; else if($3=="gene") print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "gene_id" " \"" $9 "\""; else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' annotated.gtf | sed "s/\"t/t/g" | sed "s/\"\"/\"/g" | sed "s/\;\"/;/g" > final_annotated.gtf
#awk -F"\t" '{if($3=="transcript") print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "transcript_id" " \"" $9 "\"; Name \"pfam_id>" $11 "_UniProt_id>" $14 "," $15 "," $16 "_tmhmm>" $18 "," $19 "_" $21 "\""; else if($3=="gene") print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "gene_id" " \"" $9 "\""; else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' annotated.gff | sed "s/\"t/t/g" | sed "s/\"\"/\"/g" | sed "s/\;\"/;/g" > final_annotated.gff
printf '##gff-version 3\n##Index-subfeatures 1' > final_annotated_v3.gff3 
awk -F"\t" '{split($9,a,".");split($9,b,";");split(b[1],c,"\"");if($3=="transcript") print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "ID=" $9 ";" "Parent=" a[1] ";Note=" $11 "-" $12 "," $19 "," $20 "," $21 "," $15 "," $17 ; else if($3=="gene") print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "ID=" $9 ; else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "Parent=" c[4]}' annotated.gff >> final_annotated_v3.gff3
# delete "Topology=" bcoz it contains seconds "=" sign which is not allowed
# sort final_annotated.gff
#bedtools sort -i final_annotated.gff > final_annotated_sorted.gff
bedtools sort -header -i final_annotated_v3.gff3 > final_annotated_v3_sorted.gff3

#remove all the genes with no pfam domain
#first generate list of such genes and respective transcripts
#grep 'pfam_id>_' final_annotated_sorted.gff | awk -F"\t" '{split($9,a,";");split(a[1],b,"\"");split(b[2],c,"."); print "gene_id \""c[1]"\"" "\n" "transcript_id \""b[2]"\""}' | sort > no_pfam_id_genes.lst
#list genes with no uniprot ids
#grep 'UniProt_id>,' final_annotated_sorted.gff | awk -F"\t" '{split($9,a,";");split(a[1],b,"\"");split(b[2],c,"."); print "gene_id \""c[1]"\"" "\n" "transcript_id \""b[2]"\""}' | sort > no_uniprot_id_genes.lst 
#keep genes which are common in both lists
#comm -12 no_pfam_id_genes.lst no_uniprot_id_genes.lst | sort > gene_noPfam_noUniprot.lst
#use this list to remove all matching lines
#grep -Fv -f gene_noPfam_noUniprot.lst final_annotated_sorted.gff > final_annotated_withPfamDomain_sorted.gff3

## second approach : using awk
awk -F"\t" '{split($9,a,";");split(a[1],a1,"=");split(a[2],a2,"=");split(a[3],a3,"=");split(a3[2],b,","); if((b[1]=="") && (b[2]=="") && ($3=="transcript")) print a1[2]}' final_annotated_v3_sorted.gff3 > transcripts_withnopfam_uniprotdomain.lst
grep -Fwv -f transcripts_withnopfam_uniprotdomain.lst final_annotated_v3_sorted.gff3 > final_annotated_v3_withPfamUniprotDomain_partialDone_sorted.gff3
awk -F"\t" '{split($9,a,";");split(a[1],a1,"=");split(a[2],a2,"=");if($3=="transcript") print a2[2]}' final_annotated_v3_withPfamUniprotDomain_partialDone_sorted.gff3 | sort > transcripts_withPfam_Uniprotdomain_geneIDS.lst
awk -F"\t" '{split($9,a,";");split(a[1],a1,"=");if($3=="gene") print a1[2]}' final_annotated_v3_withPfamUniprotDomain_partialDone_sorted.gff3 | sort > all_geneIDs.lst
comm -23 all_geneIDs.lst transcripts_withPfam_Uniprotdomain_geneIDS.lst > genes_presentonlywith_one_transcript_and_no_pfam_and_uniprot_domain.lst
grep -Fwv -f genes_presentonlywith_one_transcript_and_no_pfam_and_uniprot_domain.lst final_annotated_v3_withPfamUniprotDomain_partialDone_sorted.gff3 > final_annotated_v3_withPfamUniprotDomain_Done_sorted.gff3
sed -i 's/Topology=/Topology_/g' final_annotated_v3_withPfamUniprotDomain_Done_sorted.gff3
sed -i 's/initial/exon/g' final_annotated_v3_withPfamUniprotDomain_Done_sorted.gff3
sed -i 's/terminal/exon/g' final_annotated_v3_withPfamUniprotDomain_Done_sorted.gff3
sed -i 's/internal/exon/g' final_annotated_v3_withPfamUniprotDomain_Done_sorted.gff3
sed -i 's/single/exon/g' final_annotated_v3_withPfamUniprotDomain_Done_sorted.gff3

grep -vP '\tCDS\t' final_annotated_v3_withPfamUniprotDomain_Done_sorted.gff3 > express_v3_annotation.gff3

gffread express_v3_annotation.gff3 -O -F -o- | awk -F"\t" '{split($9,a,";"); split(a[1],a1,"=");split(a1[2],a12,".");if(($3=="mRNA") && ($9 !~ /geneID/)) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6  "\t" $7 "\t" $8 "\t" "geneID="a12[1]";"$9 "\t1"; else if(($3=="mRNA") && ($9 !~ /Parent/)) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6  "\t" $7 "\t" $8 "\t" "Parent="a12[1]";"$9 "\t2"; else if(($3=="mRNA") && ($9 !~ /Parent/) && ($9 !~ /geneID/)) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6  "\t" $7 "\t" $8 "\t" "geneID="a12[1]";Parent="a12[1]";"$9 "\t3"; else print $0 "\t4"}' > express_v3_annotation_gffread_awk.gff3



#########################
#filter gff3
####single exon genes with no uniprot match
grep -P '\texon\t' ../express_v3_annotation_gffread_awk.gff3 | awk '{print $9}' | uniq -c | awk '{split($2,a,"Parent="); if($1==1) print a[2]}' | grep -f - ../test_v3.gff3 | awk '{split($9,a,"Note="); split(a[2],a2,","); if((a2[2]=="")&&($3=="mRNA")) print $0}' | wc -l 



####for new blast 

printf '##gff-version 3\n##Index-subfeatures 1' > final_annotated_v4.gff3 
awk -F"\t" '{split($9,a,".");split($9,b,";");split(b[1],c,"\"");if($3=="transcript") print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "ID=" $9 ";" "Parent=" a[1] ";Note=" $11 "," $19 "," $20 "," $21 "," $22 "," $14 "," $16 ; else if($3=="gene") print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "ID=" $9 ; else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "Parent=" c[4]}' annotated.gff >> final_annotated_v4.gff3
bedtools sort -header -i final_annotated_v4.gff3 > final_annotated_v4_sorted.gff3
sed 's/,PF/-PF/g' final_annotated_v4_sorted.gff3 > final_annotated_v4.1_sorted_.gff3 #to fix comma between two pfam domains
awk -F"\t" '{split($9,a,";");split(a[1],a1,"=");split(a[2],a2,"=");split(a[3],a3,"=");split(a3[2],b,","); if((b[1]=="") && (b[2]=="") && ($3=="transcript")) print a1[2]}' final_annotated_v3_sorted.gff3 > transcripts_withnopfam_uniprotdomain.lst

# here we are getting new gene names, so I decided to include complete uniprot blast results within the same annotation out as version v4.2.1
#some changes were needed in the python script in order to process above data. This time I have included e-values too.
python3 /media/pflanz/Hs1-2/Scripts/combine_annotation_gff_toNewGFF/combine_all_toGff.py #gave us final_annotated_v4.2_sorted.gff3

awk -F"\t" '{if($3=="transcript") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9","$11","$12","$13","$14 ; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9 }' final_annotated_v4.2_sorted.gff3 > final_annotated_v4.2.1_sorted.gff3



### count genes with filter of evalue
awk -F"\t" '{split($9,a,"Note=");split(a[2],a2,","); if(($3=="transcript") && ((a2[5]<0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001))) print $0}' final_annotated_v4.2.1_sorted.gff3 | awk -F"\t" '{split($9,a,"ID=");split(a[2],a2,";");split(a2[1],a21,".");print a21[1]}' | sort | uniq | wc -l
## filter evalue and search genes
awk -F"\t" '{split($9,a,"Note=");split(a[2],a2,","); if(($3=="transcript") && ((a2[5]<0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001))) print $0}' final_annotated_v4.2.1_sorted.gff3 | awk -F"\t" '{split($9,a,"ID=");split(a[2],a2,";");split(a2[1],a21,".");print $0 "\t" a21[1]}' | grep ',SPL9_'
## find single exon genes
grep -P '\texon\t' final_annotated_v4.2.2_sorted.gff3 | awk '{split($9,a,"Parent=");split(a[2],a2,".");print a2[1]}' | sort | uniq -c | awk '{if($1==1) print $2}' > single_Exon_geneIDs_final_annotated_v4.2.2_sorted.gff3
## single exon genes with size smaller with ARATH match
grep -w -f single_Exon_geneIDs_final_annotated_v4.2.2_sorted.gff3 final_annotated_v4.2.2_sorted.gff3 | awk '{split($9,a,"Note=");split(a[2],a2,",");if($3=="transcript") print $0}' | awk '{split($9,a,"Note=");split(a[2],a2,",");if(a2[8] ~ '/_ARATH/')print $0 }' | awk '{split($9,a,"ID=");split(a[2],a2,";"); split(a2[1],a21,"."); print a21[1]}' > single_exon_geneIDs_with_ARATH_final_annotated_v4.2.2_sorted.gff3

### we want to blast against bra ol, ra, na, sinapus alba
# make new blast db with just these four
cat uniprot-prot-Brnapus.fasta uniprot-prot-BrOl.fasta uniprot-prot-BrRapa.fasta uniprot-Sinapis+alba.fasta > Rapa_Ol_Na_Sal_uniprot.fasta

########################################################################
#####################  final  ##########################################
########################################################################

#####new blast -2 saved in a folder named new_blast_2 in side annotation folder
#make blast db
~/ncbi-blast-2.7.1+/bin/makeblastdb -in Rapa_Ol_Na_Sal_uniprot.fasta -dbtype prot
#just for ARTH
~/ncbi-blast-2.7.1+/bin/makeblastdb -in ../../../../databases/uniprot-filtered-organism-Arabidopsis+thaliana5B370.fasta -dbtype prot
#blast our prot seq against above generated database
~/ncbi-blast-2.7.1+/bin/blastp -query ../../trinotate/protein_eachCDS.fasta -db ../../../../databases/Rapa_Ol_Na_Sal_uniprot.fasta -num_threads 12 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > proteinCDS_Rapa_Ol_Na_Sal_uniprot.outfmt6
# blast against uniprot ARATH only
~/ncbi-blast-2.7.1+/bin/blastp -query ../../trinotate/protein_eachCDS.fasta -db ../../../../databases/uniprot-filtered-organism-Arabidopsis+thaliana5B370.fasta -num_threads 12 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > proteinCDS_At_uniprot.outfmt6
#filter blast results
#keep only important columns from outfmt6.
#get one pfam domain for one region using below perl scipt
~/softwares/PfamScan/pfam_scan.pl -cpu 12 -dir ~/softwares/PfamScan -fasta /media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/trinotate/protein_eachCDS.fasta -outfile /media/pflanz/Avneesh/Brassicanapus/gatk/second_Attempt/annotation/pfam_scan/pfam_CDSprotein_result.out
# for v4.3
python3 /media/pflanz/Hs1-2/Scripts/combine_annotation_gff_toNewGFF/combine_all_toGff.py
##clean above generated output file
printf '##gff-version 3\n##Index-subfeatures 1\n' > final_annotated_v4.3.1.gff3
awk -F"\t" '{split($9,a,".");split($9,b,";");split(b[1],c,"\"");if($3=="transcript") print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "ID=" $9 ";" "Parent=" a[1] ";Note=" $11 "," $14 "," $15 "," $16 "," $19 "," $20 "," $21 "," $26; else if($3=="gene") print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "ID=" $9 ; else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "Parent=" c[4]}' final_annotated_v4.3.gff3 >> final_annotated_v4.3.1.gff3

#count all unique number of genes after filter for evalue blast of 10e-50 and min size of 200bp
#awk -F"\t" '{split($9,a,"Note=");split(a[2],a2,",");split($9,b,"Parent=");split(b[2],b2,";"); if(($3=="transcript")&&(a2[6]!="")&&(a2[9]<0.000000000000000000000000000000000000000000000000000000000001)&&(($5-$4)>200)) print b2[1]}' final_annotated_v4.3.1.gff3 | sort | uniq | wc -l #gives us 103091
awk -F"\t" '{split($9,a,"Note=");split(a[2],a2,",");split($9,b,"Parent=");split(b[2],b2,";"); if(($3=="transcript")&&(a2[6]!="")&&(a2[9]<0.00000000000000000000000000000000000000000000000001)) print b2[1]}' final_annotated_v4.3.1.gff3 > final_annotated_v4.3.1._filtered_geneIDs.lst
grep -Fw -f final_annotated_v4.3.1._filtered_geneIDs.lst final_annotated_v4.3.1.gff3 > final_annotated_v4.3.2.gff3
#sed 's/Topology=/Topology_/g' final_annotated_v4.3.3.gff3
sed 's/initial/exon/g' final_annotated_v4.3.2.gff3 > final_annotated_v4.3.3.gff3
sed -i 's/terminal/exon/g' final_annotated_v4.3.3.gff3
sed -i 's/internal/exon/g' final_annotated_v4.3.3.gff3
sed -i 's/single/exon/g' final_annotated_v4.3.3.gff3
grep -vP '\tCDS\t' final_annotated_v4.3.3.gff3 > final_annotated_v4.3.4.gff3
bedtools sort -header -i final_annotated_v4.3.4.gff3 > final_annotated_v4.3.4_sorted.gff3
~/softwares/gffread/gffread/gffread final_annotated_v4.3.4_sorted.gff3 -O -F -o- > final_annotated_v4.3.5_sorted.gff3
###from now on things will be run on server
# step 1: STAR mapping using 
# ste 2: sort, convert to bam
# step 3: cufflinks
# step 4: cuffdiff

### make promoter borders ####
grep -v '^gene*' ../gene_exp_analysis/v1_4_steep_up_vlTOvh.out | awk '{print $1}' | sort | uniq | grep -Fw -f -  ../annotation/final_annotated_v4.3.5_sorted.gff3 | awk '{if(($3=="mRNA") && ($7=="+")) print $1 "\t" $2 "\t" "Promoter" "\t" $4-3000 "\t" $4 "\t" $6 "\t" $7 "\t" $8 "\t" $9; else if(($3=="mRNA") && ($7=="-")) print  $1 "\t" $2 "\t" "Promoter" "\t" $5 "\t" $5+3000 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' > v1_4_steep_up_vlTOvh.out.promoter

##### get sequence for region of interest #####
bedtools getfasta -s -fi ../Express_v1.fa -bed v1_3_steep_up_vlTOvh_l.out.promoter -fo v1_3_steep_up_vlTOvh_l.out.promoter_strandness.fasta

## use meme in parallel mode to find motifs
~/softwares/meme-5.0.4/src/parallel/meme v1_4_steep_up_vlTOvh.out.promoter_stradness.fasta -o ./meme_strandess -p 12 -minw 3 -maxw 10 -dna


