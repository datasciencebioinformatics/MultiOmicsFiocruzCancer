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



# Histogram overlaid with kernel density curve
plot_enropy_stage_I   <-ggplot(df_enropy_stage_I, aes(x=entropy)) +  geom_histogram(aes(y=..density..), binwidth=.20, colour="black", fill="white") +  geom_density(alpha=.2, fill="#FF6666") + geom_segment(aes(x=entropy_stage_I, y=200, xend=entropy_stage_I, yend=0), arrow = arrow(length=unit(0.5, 'cm')))
plot_enropy_stage_II  <-ggplot(df_enropy_stage_II, aes(x=entropy)) +  geom_histogram(aes(y=..density..), binwidth=.20, colour="black", fill="white") +  geom_density(alpha=.2, fill="#FF6666")  + geom_segment(aes(x=entropy_stage_I, y=200, xend=entropy_stage_I, yend=0), arrow = arrow(length=unit(0.5, 'cm')))
plot_enropy_stage_III  <-ggplot(df_enropy_stage_III, aes(x=entropy)) +  geom_histogram(aes(y=..density..), binwidth=.20, colour="black", fill="white") +  geom_density(alpha=.2, fill="#FF6666")  + geom_segment(aes(x=entropy_stage_I, y=200, xend=entropy_stage_I, yend=0), arrow = arrow(length=unit(0.5, 'cm')))
