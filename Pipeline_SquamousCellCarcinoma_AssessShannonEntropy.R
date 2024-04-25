####################################################################################################################
unstranded_data_file                <- "/home/felipe/Documentos/LungPortal/samples/unstranded.tsv"                 #
unstranded_data                     <-read.table(file = unstranded_data_file, sep = '\t', header = TRUE,fill=TRUE) #
####################################################################################################################
# Take stageI_list_of_genes
list_of_genes<-unique(merge_interactome_gene_symbol[,c("gene_id","PPI")])

# File path
file_genes_Stage_I   <-   "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_I.tsv"
file_genes_Stage_II   <-  "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_II.tsv"
file_genes_Stage_III   <- "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_III.tsv"

# Gene table
genes_Stage_I       <-read.table(file = file_genes_Stage_I, sep = '\t', header = TRUE,fill=TRUE)         #
genes_Stage_II      <-read.table(file = file_genes_Stage_II, sep = '\t', header = TRUE,fill=TRUE)#
genes_Stage_III     <-read.table(file = file_genes_Stage_III, sep = '\t', header = TRUE,fill=TRUE)

genes_Stage_I<-list_of_genes[list_of_genes$gene %in% genes_Stage_I$gene,]
genes_Stage_II<-list_of_genes[list_of_genes$gene %in% genes_Stage_II$gene,]
genes_Stage_III<-list_of_genes[list_of_genes$gene %in% genes_Stage_III$gene,]

entropy_stage_I  <-round(Entropy(genes_Stage_I$PPI, base = 2),3)
entropy_stage_II <-round(Entropy(genes_Stage_II$PPI, base = 2),3)
entropy_stage_III<-round(Entropy(genes_Stage_III$PPI, base = 2),3)

entropy_bootstrapping_stage_I     <-c()
entropy_bootstrapping_stage_II    <-c()
entropy_bootstrapping_stage_III   <-c()

# Repeat 1000 times
for (bootstrapping in 1:1000)
{
  random_genes_Stage_I<-sample(rownames(unstranded_data), dim(genes_Stage_I)[1], replace = FALSE, prob = NULL)
  random_genes_Stage_II<-sample(rownames(unstranded_data), dim(genes_Stage_II)[1], replace = FALSE, prob = NULL)
  random_genes_Stage_III<-sample(rownames(unstranded_data), dim(genes_Stage_III)[1], replace = FALSE, prob = NULL)
  
  # Take random genes
  random_genes_Stage_I    <-list_of_genes[list_of_genes$gene %in% random_genes_Stage_I, ]
  random_genes_Stage_II   <-list_of_genes[list_of_genes$gene %in% random_genes_Stage_II, ]
  random_genes_Stage_III  <-list_of_genes[list_of_genes$gene %in% random_genes_Stage_III, ]
  
  entropy_bootstrapping_stage_I<-c(entropy_bootstrapping_stage_I,round(Entropy(random_genes_Stage_I$PPI, base = 2),3))
  entropy_bootstrapping_stage_II<-c(entropy_bootstrapping_stage_II,round(Entropy(random_genes_Stage_II$PPI, base = 2),3))
  entropy_bootstrapping_stage_III<-c(entropy_bootstrapping_stage_III,round(Entropy(random_genes_Stage_III$PPI, base = 2),3))
}
# Save stages
df_enropy_stage_I  <-data.frame(1:1000,entropy=entropy_bootstrapping_stage_I,stage="Stage I")
df_enropy_stage_II <-data.frame(1:1000,entropy=entropy_bootstrapping_stage_II,stage="Stage II")
df_enropy_stage_III<-data.frame(1:1000,entropy=entropy_bootstrapping_stage_III,stage="Stage III")

p_value_stage_I<-paste("p.value: ",sum(df_enropy_stage_I$entropy>=entropy_stage_I)/1000)
p_value_stage_II<-paste("p.value: ",sum(df_enropy_stage_II$entropy>=entropy_stage_II)/1000)
p_value_stage_III<-paste("p.value: ",sum(df_enropy_stage_III$entropy>=entropy_stage_III)/1000)

# Histogram overlaid with kernel density curve
plot_enropy_stage_I     <-ggplot(df_enropy_stage_I, aes(x=entropy))  + geom_histogram() +  geom_segment(aes(x=entropy_stage_I, y=200, xend=entropy_stage_I, yend=0), arrow = arrow(length=unit(0.5, 'cm'))) + ggtitle(paste("Entropy 1000x stage I\n",p_value_stage_I,sep=""))       + theme_bw()
plot_enropy_stage_II    <-ggplot(df_enropy_stage_II, aes(x=entropy))  + geom_histogram() +  geom_segment(aes(x=entropy_stage_II, y=200, xend=entropy_stage_II, yend=0), arrow = arrow(length=unit(0.5, 'cm'))) + ggtitle(paste("Entropy 1000x stage II\n",p_value_stage_II,sep=""))    + theme_bw()
plot_enropy_stage_III   <-ggplot(df_enropy_stage_III, aes(x=entropy))  + geom_histogram() +  geom_segment(aes(x=entropy_stage_III, y=200, xend=entropy_stage_III, yend=0), arrow = arrow(length=unit(0.5, 'cm'))) + ggtitle(paste("Entropy 1000x stage I\n",p_value_stage_III,sep="")) + theme_bw()

# FindClusters_resolution
png(filename=paste(output_dir,"Entropy_","all_.png",sep=""), width = 24, height = 16, res=600, units = "cm")
	print(grid.arrange(plot_enropy_stage_I,plot_enropy_stage_II,plot_enropy_stage_III, ncol=3))
dev.off()
####################################################################################################################
n_genes<-ceiling((dim(genes_Stage_I)[1]+dim(genes_Stage_II)[1]+dim(genes_Stage_I)[1])/3)                           #
n_genes=100
####################################################################################################################
library("ggpubr")

entropy_bootstrapping_stage_I     <-c()
entropy_bootstrapping_stage_II    <-c()
entropy_bootstrapping_stage_III   <-c()

# Repeat 1000 times
for (bootstrapping in 1:1000)
{
  random_genes_Stage_I<-sample(rownames(unstranded_data),   n_genes, replace = FALSE, prob = NULL)
  random_genes_Stage_II<-sample(rownames(unstranded_data),  n_genes, replace = FALSE, prob = NULL)
  random_genes_Stage_III<-sample(rownames(unstranded_data), n_genes, replace = FALSE, prob = NULL)
  
  # Take random genes
  random_genes_Stage_I    <-list_of_genes[list_of_genes$gene %in% random_genes_Stage_I, ]
  random_genes_Stage_II   <-list_of_genes[list_of_genes$gene %in% random_genes_Stage_II, ]
  random_genes_Stage_III  <-list_of_genes[list_of_genes$gene %in% random_genes_Stage_III, ]
  
  entropy_bootstrapping_stage_I<-c(entropy_bootstrapping_stage_I,round(Entropy(random_genes_Stage_I$PPI, base = 2),3))
  entropy_bootstrapping_stage_II<-c(entropy_bootstrapping_stage_II,round(Entropy(random_genes_Stage_II$PPI, base = 2),3))
  entropy_bootstrapping_stage_III<-c(entropy_bootstrapping_stage_III,round(Entropy(random_genes_Stage_III$PPI, base = 2),3))
}
# Entropy bootstrapping
entropy_bootstrapping<-rbind(data.frame(Stage="Stage I",entropy=entropy_bootstrapping_stage_I),data.frame(Stage="Stage II",entropy=entropy_bootstrapping_stage_II),data.frame(Stage="Stage III",entropy=entropy_bootstrapping_stage_III))

# My comparissons
my_comparisons <- list( c("Stage I", "Stage II"), c("Stage I", "Stage III"), c("Stage II", "Stage III") )

# Create plot
entropy_distribution_plot<-ggplot(entropy_bootstrapping, aes(x=entropy,fill=Stage))  + geom_density(alpha=.3)  + ggtitle(paste("Entropy densities 1000x all stages",sep="")) + theme_bw() +  scale_fill_manual(values = c("black", "orange", "purple"))
entropy_boxplot_plot<-ggboxplot(entropy_bootstrapping, x = "Stage", y = "entropy", color = "Stage") + stat_compare_means(comparisons = my_comparisons) + theme_bw() +  scale_colour_manual(values = c("black", "orange", "purple"))+ ggtitle(paste("Entropy boxplot 1000x all stages",sep=""))

# FindClusters_resolution
png(filename=paste(output_dir,"Entropy_","distribution_boxplot.png",sep=""), width = 28, height = 16, res=600, units = "cm")
	print(grid.arrange(entropy_distribution_plot,entropy_boxplot_plot, ncol=2))
dev.off()
