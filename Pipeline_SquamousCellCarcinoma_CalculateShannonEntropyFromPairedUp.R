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
file_genes_Stage_I   <-   paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_I.tsv",sep="")
file_genes_Stage_II   <-  paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_II.tsv",sep="")
file_genes_Stage_III   <- paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_III.tsv",sep="")

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
########################################################################################################################################
df_stageI_connectivity   <-unique(data.frame(Conectivity=table(c(interactome_stage_I$Gene1,interactome_stage_I$Gene2))))
df_stageII_connectivity  <-unique(data.frame(Conectivity=table(c(interactome_stage_II$Gene1,interactome_stage_II$Gene2))))
df_stageIII_connectivity <-unique(data.frame(Conectivity=table(c(interactome_stage_III$Gene1,interactome_stage_III$Gene2))))
########################################################################################################################################
colnames(df_stageI_connectivity)<-c("Gene","Conectivity")
colnames(df_stageII_connectivity)<-c("Gene","Conectivity")
colnames(df_stageIII_connectivity)<-c("Gene","Conectivity")
########################################################################################################################################

# Table for the calculation of entropy
df_entropy_calulation_I   <-data.frame(table(df_stageI_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
df_entropy_calulation_II  <-data.frame(table(df_stageII_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
df_entropy_calulation_III <-data.frame(table(df_stageIII_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)

# Rename colnames
colnames(df_entropy_calulation_I)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
colnames(df_entropy_calulation_II)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
colnames(df_entropy_calulation_III)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")

# Calculate p(k)
df_entropy_calulation_I$p_k<-df_entropy_calulation_I$count/sum(df_entropy_calulation_I$count)
df_entropy_calulation_II$p_k<-df_entropy_calulation_II$count/sum(df_entropy_calulation_II$count)
df_entropy_calulation_III$p_k<-df_entropy_calulation_III$count/sum(df_entropy_calulation_III$count)

# Calculate log2(p(k))
df_entropy_calulation_I$log2_pk<-log(df_entropy_calulation_I$p_k,2)
df_entropy_calulation_II$log2_pk<-log(df_entropy_calulation_II$p_k,2)
df_entropy_calulation_III$log2_pk<-log(df_entropy_calulation_III$p_k,2)

# Calculate p(k)*log2(p(k))
df_entropy_calulation_I$p_k_mult_log2_pk<-df_entropy_calulation_I$p_k*df_entropy_calulation_I$log2_pk
df_entropy_calulation_II$p_k_mult_log2_pk<-df_entropy_calulation_II$p_k*df_entropy_calulation_II$log2_pk
df_entropy_calulation_III$p_k_mult_log2_pk<-df_entropy_calulation_III$p_k*df_entropy_calulation_III$log2_pk

# Caclulate entropy value
Entropy_stage_I_value_Carels  <-abs(sum(df_entropy_calulation_I$p_k_mult_log2_pk))
Entropy_stage_II_value_Carels <-abs(sum(df_entropy_calulation_II$p_k_mult_log2_pk))
Entropy_stage_III_value_Carels<-abs(sum(df_entropy_calulation_III$p_k_mult_log2_pk))
########################################################################################################################################
# Save TSV file with genes from Stage1
write_tsv(df_stageI_connectivity, paste(output_dir,"df_stageI_connectivity_I_interactome",".tsv",sep=""))
write_tsv(df_stageII_connectivity, paste(output_dir,"df_stageII_connectivity_II_interactome",".tsv",sep=""))
write_tsv(df_stageIII_connectivity, paste(output_dir,"df_stageIII_connectivity_III_interactome",".tsv",sep=""))

# Save TSV file with genes from Stage1
write_tsv(interactome_data_stage_I, paste(output_dir,"df_stageI_interactome_interactome",".tsv",sep=""))
write_tsv(interactome_data_stage_II, paste(output_dir,"df_stageII_interactome_interactome",".tsv",sep=""))
write_tsv(interactome_data_stage_III, paste(output_dir,"df_stageIII_interactome_interactome",".tsv",sep=""))
########################################################################################################################################
cat(print(paste("\nNº of vertex/Nº/Entropy of edges, sub-interactome network for Stage I: ",  paste(length(df_stageI_connectivity$Gene),dim(unique(interactome_stage_I))[1],round(Entropy_stage_I_value_Carels,4),sep="/"),sep="")),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
cat(print(paste("\nNº of vertex/Nº/Entropy of edges, sub-interactome network for Stage II: ", paste(length(df_stageII_connectivity$Gene),dim(unique(interactome_stage_II))[1],round(Entropy_stage_II_value_Carels,4),sep="/"),sep="")),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
cat(print(paste("\nNº of vertex/Nº/Entropy of edges, sub-interactome network for Stage III: ",paste(length(df_stageIII_connectivity$Gene),dim(unique(interactome_stage_III))[1],round(Entropy_stage_III_value_Carels,4),sep="/"),sep="")),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
