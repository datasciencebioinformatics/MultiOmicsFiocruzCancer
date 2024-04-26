library(DescTools)

# Interactome
interactome_file<-"/home/felipe/Documentos/LungPortal/H-I-05.tsv"

# Gene table
interactome_data <-read.table(file = interactome_file, sep = '\t', header = FALSE,fill=TRUE)         #

# Rename collumns 
colnames(interactome_data)<-c("Gene1","Gene2")
#######################################################################################################################################
# Take stageI_list_of_genes
list_of_genes<-unique(merge_interactome_gene_symbol[,c("gene_id","PPI")])

# File path
file_genes_Stage_I   <-   "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_I.tsv"
file_genes_Stage_II   <-  "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_II.tsv"
file_genes_Stage_III   <- "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_III.tsv"

# Gene table
genes_Stage_I       <-read.table(file = file_genes_Stage_I, sep = '\t', header = TRUE,fill=TRUE)         #
genes_Stage_II      <-read.table(file = file_genes_Stage_II, sep = '\t', header = TRUE,fill=TRUE)#
genes_Stage_III     <-read.table(file = file_genes_Stage_III, sep = '\t', header = TRUE,fill=TRUE)                 #
########################################################################################################################################
# Storte genes stage I, II and III
genes_id_vector_stage_I<-c()
genes_id_vector_stage_II<-c()
genes_id_vector_stage_III<-c()

# For each gene in stage I
for (gene_id in genes_Stage_I$gene)
{
  # Store vector
  genes_id_vector_stage_I<-c(genes_id_vector_stage_I,strsplit(gene_id, split = "\\.")[[1]][1])
}
# For each gene in stage II
for (gene_id in genes_Stage_II$gene)
{
  # Store vector
  genes_id_vector_stage_II<-c(genes_id_vector_stage_II,strsplit(gene_id, split = "\\.")[[1]][1])
}
# For each gene in stage III
for (gene_id in genes_Stage_III$gene)
{
  # Store vector
  genes_id_vector_stage_III<-c(genes_id_vector_stage_III,strsplit(gene_id, split = "\\.")[[1]][1])
}
########################################################################################################################################
# Store genes
genes<-unique(c(interactome_data$Gene1,interactome_data$Gene2))

# Save genes that are in the interactome
genes_interactome_stage_I  <-genes_id_vector_stage_I[genes_id_vector_stage_I %in% genes]
genes_interactome_stage_II <-genes_id_vector_stage_II[genes_id_vector_stage_II %in% genes]
genes_interactome_stage_III<-genes_id_vector_stage_III[genes_id_vector_stage_III %in% genes]

# Gene pairs stage I
gene_starge_I<- unique(interactome_data[c(which(interactome_data$Gene1 %in% genes_interactome_stage_I),
which(interactome_data$Gene2 %in% genes_interactome_stage_I)),])

# Gene pairs stage II
gene_starge_II<- unique(interactome_data[c(which(interactome_data$Gene2 %in% genes_interactome_stage_II),
which(interactome_data$Gene2 %in% genes_interactome_stage_II)),])

# Gene pairs stage III
gene_starge_III<- unique(interactome_data[c(which(interactome_data$Gene2 %in% genes_interactome_stage_III),
which(interactome_data$Gene2 %in% genes_interactome_stage_III)),])

df_stageI_connectivity   <-data.frame(Conectivity=table(c(gene_starge_I$Gene1,gene_starge_I$Gene2)))
df_stageII_connectivity  <-data.frame(Conectivity=table(c(gene_starge_II$Gene1,gene_starge_I$Gene2)))
df_stageIII_connectivity <-data.frame(Conectivity=table(c(gene_starge_III$Gene1,gene_starge_I$Gene2)))

genes_Stage_I<-list_of_genes[list_of_genes$gene %in% genes_Stage_I$gene,]
genes_Stage_II<-list_of_genes[list_of_genes$gene %in% genes_Stage_II$gene,]
genes_Stage_III<-list_of_genes[list_of_genes$gene %in% genes_Stage_III$gene,]

round(Entropy(df_stageI_connectivity$Conectivity.Freq, base = 2),3)
round(Entropy(df_stageII_connectivity$Conectivity.Freq, base = 2),3)
round(Entropy(df_stageIII_connectivity$Conectivity.Freq, base = 2),3)
