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
#######################################################################################################################################
# Read input table
# Gene table
interactome_data <-read.table(file = interactome_file, sep = '\t', header = FALSE,fill=TRUE)         

# Gene EnsemblToUniprotKBconversionList_data
EnsemblToUniprotKBconversionList_data <-read.table(file = EnsemblToUniprotKBconversionList_file, sep = '\t', header = TRUE,fill=TRUE)         

# Rename collumns 
colnames(interactome_data)<-c("Gene1","Gene2")
#######################################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output/"   
#######################################################################################################################################
# File path to gene stages
file_genes_Stage_I   <-   "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_I.tsv"
file_genes_Stage_II   <-  "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_II.tsv"
file_genes_Stage_III   <- "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_III.tsv"

# Gene table
genes_Stage_I       <-read.table(file = file_genes_Stage_I, sep = '\t', header = TRUE,fill=TRUE)         
genes_Stage_II      <-read.table(file = file_genes_Stage_II, sep = '\t', header = TRUE,fill=TRUE)
genes_Stage_III     <-read.table(file = file_genes_Stage_III, sep = '\t', header = TRUE,fill=TRUE)   
########################################################################################################################################
# Store genes stage I, II and III
# Vectors to store gene ids from each stage
genes_id_vector_stage_I<-c()
genes_id_vector_stage_II<-c()
genes_id_vector_stage_III<-c()

# For each gene in stage I
for (gene_id in genes_Stage_I$gene)
{
  # Store gene id in the vector
  genes_id_vector_stage_I<-c(genes_id_vector_stage_I,strsplit(gene_id, split = "\\.")[[1]][1])
}
# For each gene in stage II
for (gene_id in genes_Stage_II$gene)
{
  # Store gene id in the vector
  genes_id_vector_stage_II<-c(genes_id_vector_stage_II,strsplit(gene_id, split = "\\.")[[1]][1])
}
# For each gene in stage III
for (gene_id in genes_Stage_III$gene)
{
  # Store gene id in the vector
  genes_id_vector_stage_III<-c(genes_id_vector_stage_III,strsplit(gene_id, split = "\\.")[[1]][1])
}
#######################################################################################################
# Filter tables to keep only the gene entries that are listed in the EnsemblToUniprotKBconversionList
interactome_data<-interactome_data[interactome_data$Gene1 %in% EnsemblToUniprotKBconversionList_data$SYMBOL,]
interactome_data<-interactome_data[interactome_data$Gene2 %in% EnsemblToUniprotKBconversionList_data$SYMBOL,]

# Create a table for id conversion gene_id and gene_symbol for the genes in the interactome data
gene1_conversion<-merge(interactome_data,EnsemblToUniprotKBconversionList_data,by.x="Gene1", by.y="SYMBOL",all.x=TRUE,all.y=FALSE)
gene_conversion<-merge(gene1_conversion,EnsemblToUniprotKBconversionList_data,by.x="Gene2", by.y="SYMBOL",all.x=TRUE,all.y=FALSE)

# Keep only the collumns of interest- 
# interactome_data : interactome with converted ids
interactome_data<-unique(gene_conversion[,3:4])

# Rename interactome_data collumns
colnames(interactome_data)<-c("Gene1","Gene2")
#######################################################################################################
# A filter to keep only genes that are positivelly regulated
genes_ids<-c()

# For each gene in stage I
for (gene_id in rownames(unstranded_data_filter))
{
  # Store gene id in the vector
  genes_ids<-c(genes_ids,strsplit(gene_id, split = "\\.")[[1]][1])
}
# Gene_ids
genes_ids<-unique(genes_ids)
#######################################################################################################      
# Take all genes from interactom
# store both, gene in pair one and gene in pair two in a same vectors
genes<-unique(c(interactome_data$Gene1,interactome_data$Gene2))

# The gene in lists of genes per stage must be filterd to keep only entris that are present in the interactome
# ~99% of genes in the selected lists are in the interactome
genes_interactome_stage_I  <-genes_id_vector_stage_I[genes_id_vector_stage_I %in% genes]
genes_interactome_stage_II <-genes_id_vector_stage_II[genes_id_vector_stage_II %in% genes]
genes_interactome_stage_III<-genes_id_vector_stage_III[genes_id_vector_stage_III %in% genes]
########################################################################################################################################
# If at least one of the genes in the pair are in the interactome
interactome_data_stage_I<-rbind(interactome_data[interactome_data$Gene1 %in% genes_interactome_stage_I,],
interactome_data[interactome_data$Gene2 %in% genes_interactome_stage_I,])

# If at least one of the genes in the pair are in the interactome
interactome_data_stage_II<-rbind(interactome_data[interactome_data$Gene1 %in% genes_interactome_stage_II,],
interactome_data[interactome_data$Gene2 %in% genes_interactome_stage_II,])

# If at least one of the genes in the pair are in the interactome
interactome_data_stage_III<-rbind(interactome_data[interactome_data$Gene1 %in% genes_interactome_stage_III,],
interactome_data[interactome_data$Gene2 %in% genes_interactome_stage_III,])
########################################################################################################################################
df_stageI_connectivity   <-unique(data.frame(Conectivity=table(c(interactome_data_stage_I$Gene1,interactome_data_stage_I$Gene2))))
df_stageII_connectivity  <-unique(data.frame(Conectivity=table(c(interactome_data_stage_II$Gene1,interactome_data_stage_II$Gene2))))
df_stageIII_connectivity <-unique(data.frame(Conectivity=table(c(interactome_data_stage_III$Gene1,interactome_data_stage_III$Gene2))))

round(Entropy(df_stageI_connectivity$Conectivity.Freq, base = 2),3)
round(Entropy(df_stageII_connectivity$Conectivity.Freq, base = 2),3)
round(Entropy(df_stageIII_connectivity$Conectivity.Freq, base = 2),3)

# Save TSV file with genes from Stage1
write_tsv(df_stageI_connectivity, paste(output_dir,"df_stageI_connectivity_I",".tsv",sep=""))
write_tsv(df_stageII_connectivity, paste(output_dir,"df_stageII_connectivity_II",".tsv",sep=""))
write_tsv(df_stageIII_connectivity, paste(output_dir,"df_stageIII_connectivity_II",".tsv",sep=""))

# Save TSV file with genes from Stage1
write_tsv(df_stageI_connectivity, paste(output_dir,"df_stageI_interactome",".tsv",sep=""))
write_tsv(df_stageII_connectivity, paste(output_dir,"df_stageII_interactome",".tsv",sep=""))
write_tsv(df_stageIII_connectivity, paste(output_dir,"df_stageIII_interactome",".tsv",sep=""))
########################################################################################################################################
