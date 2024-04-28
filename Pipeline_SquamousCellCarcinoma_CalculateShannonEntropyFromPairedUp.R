library("readr")
library(DescTools)
#######################################################################################################################################
# A script to caluclate entropy from lists of genes from each stage
#######################################################################################################################################
# Path to input files
# Interactome file
interactome_file<-"/home/felipe/Documentos/LungPortal/Full_Interactome_Flavia.txt"

# EnsemblToUniprotKBconversionList file
EnsemblToUniprotKBconversionList_file<-"/home/felipe/Documentos/LungPortal/EnsemblToUniprotKBconversionList.txt"

# File path
file_genes_Stage_I   <-   "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_I.tsv"
file_genes_Stage_II   <-  "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_II.tsv"
file_genes_Stage_III   <- "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_III.tsv"
#######################################################################################################################################
# Read input table
# Gene table
interactome_data <-read.table(file = interactome_file, sep = '\t', header = TRUE,fill=TRUE)         

# Gene EnsemblToUniprotKBconversionList_data
EnsemblToUniprotKBconversionList_data <-read.table(file = EnsemblToUniprotKBconversionList_file, sep = '\t', header = TRUE,fill=TRUE)         

# Rename collumns 
colnames(interactome_data)<-c("Gene1","Gene2")
#######################################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output/"   
#######################################################################################################################################
# Take stageI_list_of_genes
list_of_genes<-unique(merge_interactome_gene_symbol[,c("gene_id","PPI")])

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
# Filter tables
interactome_data<-interactome_data[interactome_data$Gene1 %in% EnsemblToUniprotKBconversionList_data$SYMBOL,]
interactome_data<-interactome_data[interactome_data$Gene2 %in% EnsemblToUniprotKBconversionList_data$SYMBOL,]

# Conversion table
gene1_conversion<-merge(interactome_data,EnsemblToUniprotKBconversionList_data,by.x="Gene1", by.y="SYMBOL",all.x=TRUE,all.y=FALSE)
gene_conversion<-merge(gene1_conversion,EnsemblToUniprotKBconversionList_data,by.x="Gene2", by.y="SYMBOL",all.x=TRUE,all.y=FALSE)

interactome_data<-gene_conversion[,3:4]
colnames(interactome_data)<-c("Gene1","Gene2")
      
# Store genes
genes<-unique(c(interactome_data$Gene1,interactome_data$Gene2))

# Save genes that are in the interactome
genes_interactome_stage_I  <-genes_id_vector_stage_I[genes_id_vector_stage_I %in% genes]
genes_interactome_stage_II <-genes_id_vector_stage_II[genes_id_vector_stage_II %in% genes]
genes_interactome_stage_III<-genes_id_vector_stage_III[genes_id_vector_stage_III %in% genes]

# Gene pairs stage I
gene_starge_I<- gene_starge_I[gene_starge_I$Gene1 %in% genes_interactome_stage_I,]
gene_starge_I<- gene_starge_I[gene_starge_I$Gene2 %in% genes_interactome_stage_I,]

# Gene pairs stage II
gene_starge_II<- gene_starge_II[gene_starge_II$Gene1 %in% genes_interactome_stage_II,]
gene_starge_II<- gene_starge_II[gene_starge_II$Gene2 %in% genes_interactome_stage_II,]

# Gene pairs stage III
gene_starge_III<- gene_starge_III[gene_starge_III$Gene1 %in% genes_interactome_stage_III,]
gene_starge_III<- gene_starge_III[gene_starge_III$Gene2 %in% genes_interactome_stage_III,]

df_stageI_connectivity   <-data.frame(Conectivity=table(c(gene_starge_I$Gene1,gene_starge_I$Gene2)))
df_stageII_connectivity  <-data.frame(Conectivity=table(c(gene_starge_II$Gene1,gene_starge_I$Gene2)))
df_stageIII_connectivity <-data.frame(Conectivity=table(c(gene_starge_III$Gene1,gene_starge_I$Gene2)))

round(Entropy(df_stageI_connectivity$Conectivity.Freq, base = 2),3)
round(Entropy(df_stageII_connectivity$Conectivity.Freq, base = 2),3)
round(Entropy(df_stageIII_connectivity$Conectivity.Freq, base = 2),3)

# Save TSV file with genes from Stage1
write_tsv(df_stageI_connectivity, paste(output_dir,"df_stageI_connectivity_I",".tsv",sep=""))
write_tsv(df_stageII_connectivity, paste(output_dir,"df_stageII_connectivity_II",".tsv",sep=""))
write_tsv(df_stageIII_connectivity, paste(output_dir,"df_stageIII_connectivity_II",".tsv",sep=""))
