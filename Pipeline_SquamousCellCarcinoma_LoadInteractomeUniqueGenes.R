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
# Select only tumor metadata
colDta_tumor<-colData[colData$tissue_type=="Tumor",]

# Select only tumor samples id's
tumor_samples<-colDta_tumor$patient_id

# Select samples of each stage stored in colData                                                                                             #
sample_stage_I  <-colDta_tumor[colDta_tumor$stages=="Stage I","patient_id"]                                                                  #
sample_stage_II <-colDta_tumor[colDta_tumor$stages=="Stage II","patient_id"]                                                                 #
sample_stage_III<-colDta_tumor[colDta_tumor$stages=="Stage III","patient_id"]                                                                #
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

interactome_stage_I$Stage<-"Stage I"
interactome_stage_II$Stage<-"Stage II"
interactome_stage_III$Stage<-"Stage III"
#######################################################################################################
interactome_all_stage<-rbind(interactome_stage_I,interactome_stage_II,interactome_stage_III)

# Remove edge
interactome_all_stage<-interactome_all_stage[interactome_all_stage$Gene1!="REMOVE",]
#######################################################################################################
df_correlation_net_stage_I<-data.frame(na.omit(unstranded_data_filter[genes_Stage_I$gene,sample_stage_I]))
df_correlation_net_stage_II<-data.frame(na.omit(unstranded_data_filter[genes_Stage_II$gene,sample_stage_II]))
df_correlation_net_stage_III<-data.frame(na.omit(unstranded_data_filter[genes_Stage_III$gene,sample_stage_III]))
#######################################################################################################################################
# Filter by low variability
# Incosistency of low-variability genes.
# Set threshold
upper_weight_th = 0.85

net_stage_I   <- cor(t(df_correlation_net_stage_I), method = "pearson", use = "complete.obs")
net_stage_II   <- cor(t(df_correlation_net_stage_II), method = "pearson", use = "complete.obs")
net_stage_III   <- cor(t(df_correlation_net_stage_III), method = "pearson", use = "complete.obs")

net_stage_I[lower.tri(net_stage_I)] <- 0
net_stage_II[lower.tri(net_stage_II)] <- 0
net_stage_III[lower.tri(net_stage_III)] <- 0

diag(net_stage_I)<-0
diag(net_stage_II)<-0
diag(net_stage_III)<-0

net_stage_I_correlation_network<-melt(net_stage_I)
net_stage_II_correlation_network<-melt(net_stage_II)
net_stage_III_correlation_network<-melt(net_stage_III)


net_stage_I_correlation_network<-na.omit(net_stage_I_correlation_network[abs(net_stage_I_correlation_network$value)>=upper_weight_th,])
net_stage_II_correlation_network<-na.omit(net_stage_II_correlation_network[abs(net_stage_II_correlation_network$value)>=upper_weight_th,])
net_stage_III_correlation_network<-na.omit(net_stage_III_correlation_network[abs(net_stage_III_correlation_network$value)>=upper_weight_th,])

# If there is not interaction
if (dim(net_stage_I_correlation_network)[1]==0)
{
  # Add edge to be removed
  net_stage_I_correlation_network<-data.frame(Gene1="REMOVE",Gene2="REMOVE",cor=0)
}
# If there is not interaction
if (dim(net_stage_II_correlation_network)[1]==0)
{
  # Add edge to be removed
  net_stage_II_correlation_network<-data.frame(Gene1="REMOVE",Gene2="REMOVE",cor=0)
}
# If there is not interaction
if (dim(net_stage_III_correlation_network)[1]==0)
{
  # Add edge to be removed
  net_stage_III_correlation_network<-data.frame(Gene1="REMOVE",Gene2="REMOVE",cor=0)
}

net_stage_I_correlation_network$Stage<-"Stage I"
net_stage_II_correlation_network$Stage<-"Stage II"
net_stage_III_correlation_network$Stage<-"Stage III"

colnames(net_stage_I_correlation_network)<-c("Gene1","Gene2","Correlation","Stage")
colnames(net_stage_II_correlation_network)<-c("Gene1","Gene2","Correlation","Stage")
colnames(net_stage_III_correlation_network)<-c("Gene1","Gene2","Correlation","Stage")
#######################################################################################################
net_stage_all_correlation_network<-rbind(net_stage_I_correlation_network,net_stage_II_correlation_network,net_stage_III_correlation_network)
#######################################################################################################
# Remove edge
net_stage_all_correlation_network<-net_stage_all_correlation_network[net_stage_all_correlation_network$Gene1!="REMOVE",]
#######################################################################################################
# Store genes stage I, II and III
# Vectors to store gene ids from each stage
genes_id_vector_stage_Gene1<-c()
genes_id_vector_stage_Gene2<-c()

# For each gene in stage I
for (gene_id in net_stage_all_correlation_network$Gene1)
{
  # Store gene id in the vector
  genes_id_vector_stage_Gene1<-c(genes_id_vector_stage_Gene1,strsplit(gene_id, split = "\\.")[[1]][1])
}
# For each gene in stage I
for (gene_id in net_stage_all_correlation_network$Gene2)
{
  # Store gene id in the vector
  genes_id_vector_stage_Gene2<-c(genes_id_vector_stage_Gene2,strsplit(gene_id, split = "\\.")[[1]][1])
}
net_stage_all_correlation_network$Gene1<-genes_id_vector_stage_Gene1
net_stage_all_correlation_network$Gene2<-genes_id_vector_stage_Gene2
#######################################################################################################
# Load interactome
genes_in_interactome<-unique(genes_Stage_ALL[genes_Stage_ALL$gene_id %in% intersect(genes_Stage_ALL$gene_id,c(interactome_all_stage$Gene1,interactome_all_stage$Gene2,net_stage_all_correlation_network$Gene1,net_stage_all_correlation_network$Gene2)),])

# Set rownames as gene_ids
rownames(genes_in_interactome)<-genes_in_interactome$gene_id

# Set gene symbols
interactome_all_stage$Gene1<-genes_in_interactome[interactome_all_stage$Gene1,"SYMBOL"]
interactome_all_stage$Gene2<-genes_in_interactome[interactome_all_stage$Gene2,"SYMBOL"]


net_stage_all_correlation_network$Gene1<-genes_in_interactome[net_stage_all_correlation_network$Gene1,"SYMBOL"]
net_stage_all_correlation_network$Gene2<-genes_in_interactome[net_stage_all_correlation_network$Gene2,"SYMBOL"]











