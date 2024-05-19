#######################################################################################################################################
# A script to caluclate entropy from lists of genes from each stage
# To construct the gene coexpression network, low variation genes are removed based on a threshold for the percentage of genes to be maintained (75%). Then a correlation matrix is ​​constructed using Spearman rank correlation, but only the upper diagonal is kept to avoid redundant edges. Finally, a correlation threshold of XX was used to maintain edges with significant coexpression.
#######################################################################################################################################
# Load conversion table
Table1_data<-read.table(file = "/home/felipe/Documentos/LungPortal/EnsemblToUniprotKBconversionList.txt", sep = '\t', header = TRUE,fill=TRUE)    
colnames(Table1_data)<-c("gene_id","gene_symbol")
Table1_data <- Table1_data[match(unique(Table1_data$gene_id), Table1_data$gene_id),]
rownames(Table1_data)<-Table1_data$gene_id
#######################################################################################################################################
# File path to gene stages
# Version 1
file_genes_Stage_I   <-   paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_I.tsv",sep="")
file_genes_Stage_II   <-  paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_II.tsv",sep="")
file_genes_Stage_III   <- paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_III.tsv",sep="")

# Version 2
#file_genes_Stage_I     <-   paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_Stage_","sample_stage_I",".tsv",sep="")
#file_genes_Stage_II    <-   paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_Stage_","sample_stage_II",".tsv",sep="")
#file_genes_Stage_III   <-   paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_Stage_","sample_stage_III",".tsv",sep="")

# Gene table
genes_Stage_I       <-read.table(file = file_genes_Stage_I, sep = '\t', header = TRUE,fill=TRUE)         
genes_Stage_II      <-read.table(file = file_genes_Stage_II, sep = '\t', header = TRUE,fill=TRUE)
genes_Stage_III     <-read.table(file = file_genes_Stage_III, sep = '\t', header = TRUE,fill=TRUE) 
#######################################################################################################################################
# Select only tumor
colDta_tumor<-colData[colData$tissue_type=="Tumor",]
colDta_tumor<-colData

# Samples of each stage stored in colData                                                                                             #
sample_stage_I  <-colDta_tumor[colDta_tumor$stages=="Stage I","patient_id"]                                                                     #
sample_stage_II <-colDta_tumor[colDta_tumor$stages=="Stage II","patient_id"]                                                                    #
sample_stage_III<-colDta_tumor[colDta_tumor$stages=="Stage III","patient_id"]                                                                   #
#######################################################################################################################################
df_correlation_net_stage_I<-data.frame(na.omit(unstranded_data_filter[genes_Stage_I$gene,]))
df_correlation_net_stage_II<-data.frame(na.omit(unstranded_data_filter[genes_Stage_II$gene,]))
df_correlation_net_stage_III<-data.frame(na.omit(unstranded_data_filter[genes_Stage_III$gene,]))
#######################################################################################################################################
# Filter by low variability
# Incosistency of low-variability genes.
# df_stage_I_filtered<-filter_low_var(t(df_correlation_net_stage_I), pct = 0.75, type = c("median"))
# df_stage_II_filtered<-filter_low_var(t(df_correlation_net_stage_II), pct = 0.75,type = c("median"))
# df_stage_III_filtered<-filter_low_var(t(df_correlation_net_stage_III), pct = 0.75, type = c("median"))

# Set threshold
upper_weight_th = threshold_cor

net_stage_I   <- cor(t(df_correlation_net_stage_I), method = "pearson", use = "complete.obs")
net_stage_II   <- cor(t(df_correlation_net_stage_II), method = "pearson", use = "complete.obs")
net_stage_III   <- cor(t(df_correlation_net_stage_III), method = "pearson", use = "complete.obs")

#net_stage_I   <- cor(df_stage_I_filtered, method = "pearson", use = "complete.obs")
#net_stage_II   <- cor(df_stage_II_filtered, method = "pearson", use = "complete.obs")
#net_stage_III   <- cor(df_stage_III_filtered, method = "pearson", use = "complete.obs")

net_stage_I[lower.tri(net_stage_I)] <- NA
net_stage_II[lower.tri(net_stage_II)] <- NA
net_stage_III[lower.tri(net_stage_III)] <- NA

net_stage_I_correlation_network<-melt(net_stage_I)
net_stage_II_correlation_network<-melt(net_stage_II)
net_stage_III_correlation_network<-melt(net_stage_III)

net_stage_I_correlation_network<-na.omit(net_stage_I_correlation_network[net_stage_I_correlation_network$value>=upper_weight_th,])
net_stage_II_correlation_network<-na.omit(net_stage_II_correlation_network[net_stage_II_correlation_network$value>=upper_weight_th,])
net_stage_III_correlation_network<-na.omit(net_stage_III_correlation_network[net_stage_III_correlation_network$value>=upper_weight_th,])
#######################################################################################################################################
# If at least one of the genes in the pair are in the interactome
interactome_data_stage_I<-net_stage_I_correlation_network[,1:2]

# If at least one of the genes in the pair are in the interactome
interactome_data_stage_II<-net_stage_II_correlation_network[,1:2]

# If at least one of the genes in the pair are in the interactome
interactome_data_stage_III<-net_stage_III_correlation_network[,1:2]

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
cat(print(paste("\nNº of vertex/Nº/Entropy of edges, co-expression network for Stage I: ",  paste(length(df_stageI_connectivity$Gene),dim(unique(interactome_data_stage_I))[1],round(Entropy_stage_I_value_Carels,4),sep="/"),sep="")),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
cat(print(paste("\nNº of vertex/Nº/Entropy of edges, co-expression network for Stage II: ", paste(length(df_stageII_connectivity$Gene),dim(unique(interactome_data_stage_II))[1],round(Entropy_stage_II_value_Carels,4),sep="/"),sep="")),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
cat(print(paste("\nNº of vertex/Nº/Entropy of edges, co-expression network for Stage III: ",paste(length(df_stageIII_connectivity$Gene),dim(unique(interactome_data_stage_III))[1],round(Entropy_stage_III_value_Carels,4),sep="/"),sep="")),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
