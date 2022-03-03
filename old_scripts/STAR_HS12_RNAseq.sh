#!/bin/sh

/home/planz/Downloads/STAR-2.5.3a/bin/Linux_x86_64/STAR \
--runThreadN 10 --genomeDir /media/planz/Avneesh/Beet_translocation/ProperDataSet/STARRNAseq/Genome \
--sjdbGTFtagExonParentTranscript /media/planz/Avneesh/Beet_translocation/ProperDataSet/Transcriptome/codingGeneFeatures.gff --sjdbOverhang 99 \
--readFilesIn /media/planz/Avneesh/Beet_translocation/ProperDataSet/Transcriptome/Reads/B0982_GAGTGG_L005_R1_001.fastq /media/planz/Avneesh/Beet_translocation/ProperDataSet/Transcriptome/Reads/B0982_GAGTGG_L005_R2_001.fastq \
--twopassMode Basic \
--outSAMstrandField intronMotif
