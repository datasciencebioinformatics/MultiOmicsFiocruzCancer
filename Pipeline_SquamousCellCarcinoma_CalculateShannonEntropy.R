library("readr")
library(DescTools)
library("GWENA")
library(data.table)
library(igraph)
#######################################################################################################################################
# A script to caluclate entropy from lists of genes from each stage
#######################################################################################################################################
# Load conversion table
Table1_data<-read.table(file = "/home/felipe/Documentos/LungPortal/EnsemblToUniprotKBconversionList.txt", sep = '\t', header = TRUE,fill=TRUE)    
colnames(Table1_data)<-c("gene_id","gene_symbol")
Table1_data <- Table1_data[match(unique(Table1_data$gene_id), Table1_data$gene_id),]
rownames(Table1_data)<-Table1_data$gene_id
#######################################################################################################################################
# File path to gene stages
file_genes_Stage_I   <-   paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_I.tsv",sep="")
file_genes_Stage_II   <-  paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_II.tsv",sep="")
file_genes_Stage_III   <- paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_III.tsv",sep="")

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
# Filter table
#df_stage_I_filtered <- filter_low_var(t(df_correlation_net_stage_I), pct = 0.75, type = "mean")
#df_stage_II_filtered <- filter_low_var(t(df_correlation_net_stage_II), pct = 0.75, type = "mean")
#df_stage_III_filtered <- filter_low_var(t(df_correlation_net_stage_III), pct = 0.75, type = "mean")

# calculate network
net_stage_I <-   build_net(df_stage_I_filtered, cor_func = "spearman", n_threads =1)
net_stage_II <-   build_net(df_stage_II_filtered, cor_func = "spearman", n_threads =1)
net_stage_III <-   build_net(df_stage_III_filtered, cor_func = "spearman", n_threads =1)

#module_stage_I  <- detect_modules(df_stage_I_filtered,  net_stage_I$network,  detailled_result = TRUE,  merge_threshold = 0.25)
#module_stage_II <- detect_modules(df_stage_II_filtered,  net_stage_II$network,  detailled_result = TRUE,  merge_threshold = 0.25)
#module_stage_II <- detect_modules(df_stage_III_filtered,  net_stage_III$network,  detailled_result = TRUE,  merge_threshold = 0.25)

sub_clusters_modules_stage_I<- get_sub_clusters(net_stage_I$network)
sub_clusters_modules_stage_II<- get_sub_clusters(net_stage_II$network)
sub_clusters_modules_stage_III<- get_sub_clusters(net_stage_III$network)

graph_stage_I <- build_graph_from_sq_mat(net_stage_I$network)
graph_stage_II <- build_graph_from_sq_mat(net_stage_II$network)
graph_stage_III <- build_graph_from_sq_mat(net_stage_III$network)
#######################################################################################################################################
# Set threshold
upper_weight_th = 0.9995

net_stage_I$network[lower.tri(net_stage_I$network)] <- NA
net_stage_II$network[lower.tri(net_stage_II$network)] <- NA
net_stage_III$network[lower.tri(net_stage_III$network)] <- NA

net_stage_I_correlation_network<-melt(net_stage_I$network)
net_stage_II_correlation_network<-melt(net_stage_II$network)
net_stage_III_correlation_network<-melt(net_stage_III$network)

net_stage_I_correlation_network<-na.omit(net_stage_I_correlation_network[net_stage_I_correlation_network$value>=upper_weight_th,])
net_stage_II_correlation_network<-na.omit(net_stage_II_correlation_network[net_stage_II_correlation_network$value>=upper_weight_th,])
net_stage_III_correlation_network<-na.omit(net_stage_III_correlation_network[net_stage_III_correlation_network$value>=upper_weight_th,])
#######################################################################################################################################
interactions_stage_I<-unique(net_stage_I_correlation_network[,c(1,2)])
interactions_stage_II<-unique(net_stage_II_correlation_network[,c(1,2)])
interactions_stage_III<-unique(net_stage_III_correlation_network[,c(1,2)])
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

round(Entropy_stage_II_value_Carels,4)
round(Entropy_stage_III_value_Carels,4)

# FindClusters_resolution
png(filename=paste(output_dir,"Network_","graph_Stage_I_.png",sep=""), width = 20, height = 20, res=600, units = "cm")
	layout_stage_I <- plot_module(graph_stage_I, upper_weight_th = upper_weight_th, groups = sub_clusters_modules_stage_I,vertex.label.cex = 0.5, node_scaling_max = 7,  legend_cex = 1,  title = paste("Network for genes of Stage I\nEntropy : ", round(Entropy_stage_I_value_Carels,4),sep=""))
dev.off()

# FindClusters_resolution
png(filename=paste(output_dir,"Network_","graph_Stage_II_.png",sep=""), width = 20, height = 20, res=600, units = "cm")
	layout_stage_II <- plot_module(graph_stage_II, upper_weight_th = upper_weight_th, groups = sub_clusters_modules_stage_II,vertex.label.cex = 0.5, node_scaling_max = 7,  legend_cex = 1,  title = paste("Network for genes of Stage II\nEntropy : ", round(Entropy_stage_II_value_Carels,4),sep=""))
dev.off()

# FindClusters_resolution
png(filename=paste(output_dir,"Network_","graph_Stage_III_.png",sep=""), width = 20, height = 20, res=600, units = "cm")
	layout_stage_III <- plot_module(graph_stage_III, upper_weight_th = upper_weight_th, groups = sub_clusters_modules_stage_III,vertex.label.cex = 0.5, node_scaling_max = 7,  legend_cex = 1,  title = paste("Network for genes of Stage III\nEntropy : ", round(Entropy_stage_III_value_Carels,4),sep=""))
dev.off()

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
