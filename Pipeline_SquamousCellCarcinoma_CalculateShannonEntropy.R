library("readr")
library(DescTools)
library("GWENA")
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

#######################################################################################################################################
df_correlation_net_stage_I<-t(data.frame(na.omit(unstranded_data_filter[genes_Stage_I$gene,sample_stage_I])))
df_correlation_net_stage_II<-t(data.frame(na.omit(unstranded_data_filter[genes_Stage_II$gene,sample_stage_II])))
df_correlation_net_stage_III<-t(data.frame(na.omit(unstranded_data_filter[genes_Stage_II$gene,sample_stage_III])))
#######################################################################################################################################
# Store genes stage I, II and III
# Vectors to store gene ids from each stage
genes_id_vector_stage_I<-c()
genes_id_vector_stage_II<-c()
genes_id_vector_stage_III<-c()

# For each gene in stage I
for (gene_id in colnames(df_correlation_net_stage_I))
{
  # Store gene id in the vector
  genes_id_vector_stage_I<-c(genes_id_vector_stage_I,strsplit(gene_id, split = "\\.")[[1]][1])
}
# For each gene in stage II
for (gene_id in colnames(df_correlation_net_stage_II))
{
  # Store gene id in the vector
  genes_id_vector_stage_II<-c(genes_id_vector_stage_II,strsplit(gene_id, split = "\\.")[[1]][1])
}
# For each gene in stage III
for (gene_id in colnames(df_correlation_net_stage_III))
{
  # Store gene id in the vector
  genes_id_vector_stage_III<-c(genes_id_vector_stage_III,strsplit(gene_id, split = "\\.")[[1]][1])
}
# Rename collumns
colnames(df_correlation_net_stage_I)<-genes_id_vector_stage_I
colnames(df_correlation_net_stage_II)<-genes_id_vector_stage_II
colnames(df_correlation_net_stage_III)<-genes_id_vector_stage_III
#######################################################################################################
# Select only tumor
colDta_tumor<-colData[colData$tissue_type=="Tumor",]

# Samples of each stage stored in colData                                                                                             #
sample_stage_I  <-colDta_tumor[colDta_tumor$stages=="Stage I","patient_id"]                                                                     #
sample_stage_II <-colDta_tumor[colDta_tumor$stages=="Stage II","patient_id"]                                                                    #
sample_stage_III<-colDta_tumor[colDta_tumor$stages=="Stage III","patient_id"]                                                                   #
#######################################################################################################################################
net_stage_I <- build_net(df_correlation_net_stage_I, cor_func = "spearman", n_threads =1)
net_stage_II <- build_net(df_correlation_net_stage_II, cor_func = "spearman", n_threads =1)
net_stage_III <- build_net(df_correlation_net_stage_III, cor_func = "spearman", n_threads =1)

# Take the correlation matrix
net_stage_I_cor<-net_stage_I$network
net_stage_II_cor<-net_stage_II$network
net_stage_III_cor<-net_stage_III$network

net_stage_I_cor <- melt(net_stage_I_cor)
net_stage_II_cor <- melt(net_stage_II_cor)  
net_stage_III_cor <- melt(net_stage_III_cor) 

# Take value of thhreshold
spearman_threshold<-0.95

dim(unique(net_stage_I_cor[net_stage_I_cor$value>=spearman_threshold,]))
dim(unique(net_stage_II_cor[net_stage_II_cor$value>=spearman_threshold,]))
dim(unique(net_stage_III_cor[net_stage_III_cor$value>=spearman_threshold,]))

interactions_stage_I<-unique(net_stage_I_cor[net_stage_I_cor$value>=spearman_threshold,c(1,2)])
interactions_stage_II<-unique(net_stage_II_cor[net_stage_II_cor$value>=spearman_threshold,c(1,2)])
interactions_stage_III<-unique(net_stage_III_cor[net_stage_III_cor$value>=spearman_threshold,c(1,2)])
#######################################################################################################################################
# If at least one of the genes in the pair are in the interactome
interactome_data_stage_I<-interactions_stage_I

# If at least one of the genes in the pair are in the interactome
interactome_data_stage_II<-interactions_stage_II

# If at least one of the genes in the pair are in the interactome
interactome_data_stage_III<-interactions_stage_III

# Rename columns
colnames(interactome_data_stage_I)<-c("Gene1","Gene2")
colnames(interactome_data_stage_II)<-c("Gene1","Gene2")
colnames(interactome_data_stage_III)<-c("Gene1","Gene2")
########################################################################################################################################
# copy interactome_data_stages
interactome_data_stage_I_clean<-interactome_data_stage_I
interactome_data_stage_II_clean<-interactome_data_stage_II
interactome_data_stage_III_clean<-interactome_data_stage_III

merge_interactome_data<-rbind(data.frame(Gene1=interactome_data_stage_I_clean$Gene1,Gene2=interactome_data_stage_I_clean$Gene2,Stage="Stage I"),
data.frame(Gene1=interactome_data_stage_II_clean$Gene1,Gene2=interactome_data_stage_II_clean$Gene2,Stage="Stage II"),
data.frame(Gene1=interactome_data_stage_III_clean$Gene1,Gene2=interactome_data_stage_III_clean$Gene2,Stage="Stage III"))

# Clean the tables
for (gene_pair_index in rownames(merge_interactome_data))
{
    # interactome_data_stage
    pair_gene_id_I <-as.vector(merge_interactome_data[gene_pair_index,"Gene1"])
    pair_gene_id_II<-as.vector(merge_interactome_data[gene_pair_index,"Gene2"])

    # Re-order gene ids
    if(pair_gene_id_II<pair_gene_id_I)
    {
      merge_interactome_data[gene_pair_index,"Gene1"]<-pair_gene_id_II
      merge_interactome_data[gene_pair_index,"Gene2"]<-pair_gene_id_I     
    }
    # Re-order gene ids
    if(pair_gene_id_II==pair_gene_id_I)
    {
      merge_interactome_data[gene_pair_index,"Gene2"]<-"REPEAT"
    }
}
# Take unique values
merge_interactome_data<-unique(merge_interactome_data)
########################################################################################################################################
interactome_data_stage_I<-merge_interactome_data[merge_interactome_data$Stage=="Stage I",]
interactome_data_stage_II<-merge_interactome_data[merge_interactome_data$Stage=="Stage II",]
interactome_data_stage_III<-merge_interactome_data[merge_interactome_data$Stage=="Stage III",]
########################################################################################################################################
interactome_data_stage_I<-unique(interactome_data_stage_I[,c("Gene1","Gene2")])
interactome_data_stage_II<-unique(interactome_data_stage_II[,c("Gene1","Gene2")])
interactome_data_stage_III<-unique(interactome_data_stage_III[,c("Gene1","Gene2")])
########################################################################################################################################
df_stageI_connectivity   <-unique(data.frame(Conectivity=table(c(interactome_data_stage_I$Gene1,interactome_data_stage_I$Gene2))))
df_stageII_connectivity  <-unique(data.frame(Conectivity=table(c(interactome_data_stage_II$Gene1,interactome_data_stage_II$Gene2))))
df_stageIII_connectivity <-unique(data.frame(Conectivity=table(c(interactome_data_stage_III$Gene1,interactome_data_stage_III$Gene2))))
########################################################################################################################################
colnames(df_stageI_connectivity)<-c("Gene","Conectivity")
colnames(df_stageII_connectivity)<-c("Gene","Conectivity")
colnames(df_stageIII_connectivity)<-c("Gene","Conectivity")
########################################################################################################################################
df_stageI_connectivity<-df_stageI_connectivity[df_stageI_connectivity$Gene!="REPEAT",]
df_stageII_connectivity<-df_stageII_connectivity[df_stageII_connectivity$Gene!="REPEAT",]
df_stageIII_connectivity<-df_stageIII_connectivity[df_stageIII_connectivity$Gene!="REPEAT",]
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
write_tsv(df_stageI_connectivity, paste(output_dir,"df_stageI_connectivity_I",".tsv",sep=""))
write_tsv(df_stageII_connectivity, paste(output_dir,"df_stageII_connectivity_II",".tsv",sep=""))
write_tsv(df_stageIII_connectivity, paste(output_dir,"df_stageIII_connectivity_II",".tsv",sep=""))

# Save TSV file with genes from Stage1
write_tsv(interactome_data_stage_I, paste(output_dir,"df_stageI_interactome",".tsv",sep=""))
write_tsv(interactome_data_stage_II, paste(output_dir,"df_stageII_interactome",".tsv",sep=""))
write_tsv(interactome_data_stage_III, paste(output_dir,"df_stageIII_interactome",".tsv",sep=""))
########################################################################################################################################
