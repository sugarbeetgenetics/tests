#import dataset
all_interesting_genes_updatedbysarah2 <- read.delim("C:/Users/suapp475/Downloads/all_interesting_genes_updated.txt", header=FALSE, comment.char="#")

#column names
cnames_update <- c("locusID", "chromosome", "start", "end", "stage_id", "FPKM", "Paralog", "gene")

#add_col_names
colnames(all_interesting_genes_updatedbysarah2) <- cnames_update

stages_order = c("q1","q2","q3","q4","q5","q6","q7","q8","q9","q10")

#plot_FPKM_vs_stage
all_interesting_genes_updatedbysarah2$gene <- as.factor(all_interesting_genes_updatedbysarah2$gene)
all_interesting_genes_updatedbysarah2$stage_id <- factor(all_interesting_genes_updatedbysarah2$stage_id, levels = stages_order)
all_interesting_genes_updatedbysarah2$FPKM <- as.numeric(all_interesting_genes_updatedbysarah2$FPKM)

library(ggplot2)
library(reshape2)
library(dplyr)

gene_names<-unique(all_interesting_genes_updatedbysarah2$gene)

for (gene in gene_names){
  print(gene)
  data_subset <- all_interesting_genes_updatedbysarah2[all_interesting_genes_updatedbysarah2$gene==gene,]
  new_data <- select(data_subset,stage_id,FPKM,Paralog)
  print(ggplot(data = new_data, aes(x=stage_id, y=FPKM, group = Paralog)) + geom_line(aes(colour=Paralog), size=1) + geom_point() + ggtitle(gene) + xlab("Stages") + ylab("FPKM"))
  show(plot)
}
