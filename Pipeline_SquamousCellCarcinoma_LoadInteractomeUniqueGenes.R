#######################################################################################################################################
# A script to caluclate entropy from lists of genes from each stage
# To construct sub-interactome networks, pairwise combinations of stage-specific genes are created and then filtered to keep edges overlapping the interctome.
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
# File path to gene stages
# File path to gene stages
# Version 1
file_genes_Stage_I    <-  paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_I",".tsv",sep="")
file_genes_Stage_II   <-  paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_II",".tsv",sep="")
file_genes_Stage_III  <-  paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_III",".tsv",sep="")

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

# Set rownames(interactome_data)
rownames(interactome_data)<-paste(interactome_data$Gene1,interactome_data$Gene2,sep="-")

# Invert data.frame
interactome_data_inv<-data.frame(Gene1=interactome_data$Gene2,Gene2=interactome_data$Gene1)

# Assert rownames
rownames(interactome_data_inv)<-paste(interactome_data_inv$Gene1,interactome_data_inv$Gene2,sep="-")

# combine interactomes
interactome_data_inv<-rbind(interactome_data,interactome_data_inv)
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

# The gene in lists of genes per stage must be filterd to keep only entris that are present in the interactome
# ~99% of genes in the selected lists are in the interactome
genes_interactome_stage_I  <-genes_id_vector_stage_I[genes_id_vector_stage_I %in% genes_ids]
genes_interactome_stage_II <-genes_id_vector_stage_II[genes_id_vector_stage_II %in% genes_ids]
genes_interactome_stage_III<-genes_id_vector_stage_III[genes_id_vector_stage_III %in% genes_ids]

# A vector with all genes of the interactome,full_interactome<-unique(c(genes_interactome_stage_I$Gene1,genes_interactome_stage_I$Gene2))
# Calculate all pairwise combinations of genes, without redunctancy
full_interactome_stage_I<- data.frame(expand.grid.unique(x = genes_interactome_stage_I, y = genes_interactome_stage_I,include.equals=FALSE))
full_interactome_stage_II<- data.frame(expand.grid.unique(x = genes_interactome_stage_II, y = genes_interactome_stage_II,include.equals=FALSE))
full_interactome_stage_III<- data.frame(expand.grid.unique(x = genes_interactome_stage_III, y = genes_interactome_stage_III,include.equals=FALSE))

# set colnames
colnames(full_interactome_stage_I)<-c("Gene1","Gene2")
colnames(full_interactome_stage_II)<-c("Gene1","Gene2")
colnames(full_interactome_stage_III)<-c("Gene1","Gene2")

rownames(full_interactome_stage_I)<-paste(full_interactome_stage_I$Gene1,full_interactome_stage_I$Gene2,sep="-")
rownames(full_interactome_stage_II)<-paste(full_interactome_stage_II$Gene1,full_interactome_stage_II$Gene2,sep="-")
rownames(full_interactome_stage_III)<-paste(full_interactome_stage_III$Gene1,full_interactome_stage_III$Gene2,sep="-")
#######################################################################################################
interactome_stage_I  <-interactome_data_inv[which(rownames(interactome_data_inv) %in% rownames(full_interactome_stage_I)),]
interactome_stage_II <-interactome_data_inv[which(rownames(interactome_data_inv) %in% rownames(full_interactome_stage_II)),]
interactome_stage_III<-interactome_data_inv[which(rownames(interactome_data_inv) %in% rownames(full_interactome_stage_III)),]

# If there is not interaction
if (dim(interactome_stage_I)[1]==0)
{
  # Add edge to be removed
  interactome_stage_I<-data.frame(Gene1="REMOVE",Gene2="REMOVE")
}
# If there is not interaction
if (dim(interactome_stage_II)[1]==0)
{
  # Add edge to be removed
  interactome_stage_II<-data.frame(Gene1="REMOVE",Gene2="REMOVE")
}
# If there is not interaction
if (dim(interactome_stage_III)[1]==0)
{
  # Add edge to be removed
  interactome_stage_III<-data.frame(Gene1="REMOVE",Gene2="REMOVE")
}
# If there is not interaction
if (dim(interactome_stage_all)[1]==0)
{
  # Add edge to be removed
  interactome_stage_all<-data.frame(Gene1="REMOVE",Gene2="REMOVE")
}

interactome_stage_I$Stage<-"Stage I"
interactome_stage_II$Stage<-"Stage II"
interactome_stage_III$Stage<-"Stage III"
interactome_stage_all$Stage<-"Stages all"
#######################################################################################################
interactome_all_stage<-rbind(interactome_stage_I,interactome_stage_II,interactome_stage_III)

# Remove edge
interactome_all_stage<-interactome_all_stage[interactome_all_stage$Gene1!="REMOVE",]
