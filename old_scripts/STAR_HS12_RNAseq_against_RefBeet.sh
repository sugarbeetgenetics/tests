#!/bin/sh
/home/planz/Downloads/STAR-2.5.3a/bin/Linux_x86_64/STAR \
--runThreadN 12 --runMode genomeGenerate --genomeDir ./RefBeet \
--genomeFastaFiles /media/planz/Avneesh/Beet_translocation/ProperDataSet/BLASTdb/RefBeet-1.2.fna; \



/home/planz/Downloads/STAR-2.5.3a/bin/Linux_x86_64/STAR \
--outReadsUnmapped Fastx \
--runThreadN 12 --genomeDir /media/planz/Avneesh/Beetxt_translocation/ProperDataSet/STARRNAseq/vs_RefBeet/RefBeet \
--sjdbGTFtagExonParentTranscript /media/planz/Avneesh/Beet_translocation/ProperDataSet/Reference/RefBeetSMRT.gff3 --sjdbOverhang 99 \
--readFilesIn /media/planz/Avneesh/Beet_translocation/ProperDataSet/Transcriptome/Reads/B0982_GAGTGG_L005_R1_001.fastq /media/planz/Avneesh/Beet_translocation/ProperDataSet/Transcriptome/Reads/B0982_GAGTGG_L005_R2_001.fastq \
--twopassMode Basic \
--outSAMstrandField intronMotif
