library("readr")
library(DescTools)
library("GWENA")
#######################################################################################################################################
# A script to caluclate entropy from lists of genes from each stage
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
# Select only tumor
colDta_tumor<-colData[colData$tissue_type=="Tumor",]

# Samples of each stage stored in colData                                                                                             #
sample_stage_I  <-colDta_tumor[colDta_tumor$stages=="Stage I","patient_id"]                                                                     #
sample_stage_II <-colDta_tumor[colDta_tumor$stages=="Stage II","patient_id"]                                                                    #
sample_stage_III<-colDta_tumor[colDta_tumor$stages=="Stage III","patient_id"]                                                                   #
#######################################################################################################################################
df_correlation_net_stage_I<-t(data.frame(na.omit(unstranded_data_filter[genes_Stage_I$gene,sample_stage_I])))
df_correlation_net_stage_II<-t(data.frame(na.omit(unstranded_data_filter[genes_Stage_II$gene,sample_stage_II])))
df_correlation_net_stage_III<-t(data.frame(na.omit(unstranded_data_filter[genes_Stage_III$gene,sample_stage_III])))

df_correlation_net_stage_I<-t(data.frame(na.omit(unstranded_data_filter[genes_Stage_I$gene,])))
df_correlation_net_stage_II<-t(data.frame(na.omit(unstranded_data_filter[genes_Stage_II$gene,])))
df_correlation_net_stage_III<-t(data.frame(na.omit(unstranded_data_filter[genes_Stage_III$gene,])))
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
net_stage_I <-   build_net(df_correlation_net_stage_I, cor_func = "spearman", n_threads =1)
net_stage_II <-  build_net(df_correlation_net_stage_II, cor_func = "spearman", n_threads =1)
net_stage_III <- build_net(df_correlation_net_stage_III, cor_func = "spearman", n_threads =1)

net_stage_I$network[lower.tri(net_stage_I$network, diag = FALSE)] <- 0
net_stage_II$network[lower.tri(net_stage_II$network, diag = FALSE)] <- 0
net_stage_III$network[lower.tri(net_stage_III$network, diag = FALSE)] <- 0

# Take the correlation matrix
net_stage_I_cor<-net_stage_I$network
net_stage_II_cor<-net_stage_II$network
net_stage_III_cor<-net_stage_III$network

net_stage_I_cor <- melt(net_stage_I_cor)
net_stage_II_cor <- melt(net_stage_II_cor)  
net_stage_III_cor <- melt(net_stage_III_cor) 

# Take value of thhreshold
spearman_threshold<-0.99

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

round(Entropy_stage_I_value_Carels,4)
round(Entropy_stage_II_value_Carels,4)
round(Entropy_stage_III_value_Carels,4)
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
g_stage_I<-graph_from_data_frame(interactome_data_stage_I, directed = TRUE, vertices = NULL)
g_stage_II<-graph_from_data_frame(interactome_data_stage_II, directed = TRUE, vertices = NULL)
g_stage_III<-graph_from_data_frame(interactome_data_stage_III, directed = TRUE, vertices = NULL)

plot(g_stage_I, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "Scale-free network model")
plot(g_stage_II, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "Scale-free network model")
plot(g_stage_III, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "Scale-free network model")
########################################################################################################################################
