############################################################################################################################################
output_folder<-"/home/felipe/Documentos/LungPortal/output/modules/"
############################################################################################################################################
# Store genes stage I, II and III
# Vectors to store gene ids from each stage
df_id_vector_stage_I<-data.frame(gene_id=c(),gene_symbol=c(), cp_gene_id=c())
df_id_vector_stage_II<-data.frame(gene_id=c(),gene_symbol=c(), cp_gene_id=c())
df_id_vector_stage_III<-data.frame(gene_id=c(),gene_symbol=c(), cp_gene_id=c())
############################################################################################################################################
# For each gene in stage I
for (cp_gene_id in rownames(df_correlation_net_stage_I))
{ 
  gene_id<-Table1_data[strsplit(cp_gene_id, split = "\\.")[[1]][1],"gene_id"]
  gene_symbol<-Table1_data[strsplit(cp_gene_id, split = "\\.")[[1]][1],"gene_symbol"]
    
  # Add row to table
  df_id_vector_stage_I<-rbind(df_id_vector_stage_I,data.frame(gene_id=gene_id,gene_symbol=gene_symbol, cp_gene_id=cp_gene_id))
}
# Remove NA lines
NA_rows<-rownames(df_id_vector_stage_I[is.na(df_id_vector_stage_I$gene_id),])

# Replace NA values for gene id values
df_id_vector_stage_I[NA_rows,"gene_id"]<-df_id_vector_stage_I[NA_rows,"cp_gene_id"]
df_id_vector_stage_I[NA_rows,"gene_symbol"]<-df_id_vector_stage_I[NA_rows,"cp_gene_id"]
############################################################################################################################################
# For each gene in stage II
for (cp_gene_id in rownames(df_correlation_net_stage_II))
{ 
  gene_id<-Table1_data[strsplit(cp_gene_id, split = "\\.")[[1]][1],"gene_id"]
  gene_symbol<-Table1_data[strsplit(cp_gene_id, split = "\\.")[[1]][1],"gene_symbol"]
    
  # Add row to table
  df_id_vector_stage_II<-rbind(df_id_vector_stage_II,data.frame(gene_id=gene_id,gene_symbol=gene_symbol, cp_gene_id=cp_gene_id))
}
# Remove NA lines
NA_rows<-rownames(df_id_vector_stage_II[is.na(df_id_vector_stage_I$gene_id),])

# Replace NA values for gene id values
df_id_vector_stage_II[NA_rows,"gene_id"]<-df_id_vector_stage_II[NA_rows,"cp_gene_id"]
df_id_vector_stage_II[NA_rows,"gene_symbol"]<-df_id_vector_stage_II[NA_rows,"cp_gene_id"]
############################################################################################################################################
# For each gene in stage II
for (cp_gene_id in rownames(df_correlation_net_stage_III))
{ 
  gene_id<-Table1_data[strsplit(cp_gene_id, split = "\\.")[[1]][1],"gene_id"]
  gene_symbol<-Table1_data[strsplit(cp_gene_id, split = "\\.")[[1]][1],"gene_symbol"]
    
  # Add row to table
  df_id_vector_stage_III<-rbind(df_id_vector_stage_III,data.frame(gene_id=gene_id,gene_symbol=gene_symbol, cp_gene_id=cp_gene_id))
}
# Remove NA lines
NA_rows<-rownames(df_id_vector_stage_III[is.na(df_id_vector_stage_I$gene_id),])

# Replace NA values for gene id values
df_id_vector_stage_III[NA_rows,"gene_id"]<-df_id_vector_stage_III[NA_rows,"cp_gene_id"]
df_id_vector_stage_III[NA_rows,"gene_symbol"]<-df_id_vector_stage_III[NA_rows,"cp_gene_id"]
############################################################################################################################################
df_id_vector_stage_I<-setDT(df_id_vector_stage_I)[, dupID := rowid(gene_symbol)]
df_id_vector_stage_II<-setDT(df_id_vector_stage_II)[, dupID := rowid(gene_symbol)]
df_id_vector_stage_III<-setDT(df_id_vector_stage_III)[, dupID := rowid(gene_symbol)]

# Create indexed values
df_id_vector_stage_I$dupID<-paste(df_id_vector_stage_I$gene_symbol,df_id_vector_stage_I$dupID-1,sep=".")
df_id_vector_stage_II$dupID<-paste(df_id_vector_stage_II$gene_symbol,df_id_vector_stage_II$dupID-1,sep=".")
df_id_vector_stage_III$dupID<-paste(df_id_vector_stage_III$gene_symbol,df_id_vector_stage_III$dupID-1,sep=".")

# remove .0 values
df_id_vector_stage_I$dupID<-gsub("\\.0", "", df_id_vector_stage_I$dupID)
df_id_vector_stage_II$dupID<-gsub("\\.0", "", df_id_vector_stage_II$dupID)
df_id_vector_stage_III$dupID<-gsub("\\.0", "", df_id_vector_stage_III$dupID)

rownames(df_id_vector_stage_I)<-df_id_vector_stage_I$cp_gene_id
rownames(df_id_vector_stage_II)<-df_id_vector_stage_II$cp_gene_id
rownames(df_id_vector_stage_III)<-df_id_vector_stage_III$cp_gene_id

df_id_vector_stage_I<-data.frame(df_id_vector_stage_I)
df_id_vector_stage_II<-data.frame(df_id_vector_stage_II)
df_id_vector_stage_III<-data.frame(df_id_vector_stage_III)

rownames(df_id_vector_stage_I)<-df_id_vector_stage_I$cp_gene_id
rownames(df_id_vector_stage_II)<-df_id_vector_stage_II$cp_gene_id
rownames(df_id_vector_stage_III)<-df_id_vector_stage_III$cp_gene_id

rownames(df_correlation_net_stage_I)<-df_id_vector_stage_I[rownames(df_correlation_net_stage_I),"dupID"]
rownames(df_correlation_net_stage_II)<-df_id_vector_stage_II[rownames(df_correlation_net_stage_II),"dupID"]
rownames(df_correlation_net_stage_III)<-df_id_vector_stage_III[rownames(df_correlation_net_stage_III),"dupID"]
############################################################################################################################################
# Filter table
df_stage_I_filtered <- filter_low_var(t(df_correlation_net_stage_I), pct = 0.75, type = "mean")
df_stage_II_filtered <- filter_low_var(t(df_correlation_net_stage_II), pct = 0.75, type = "mean")
df_stage_III_filtered <- filter_low_var(t(df_correlation_net_stage_III), pct = 0.75, type = "mean")

# calculate network
net_stage_I <-   build_net(df_stage_I_filtered, cor_func = "spearman", n_threads =1)
net_stage_II <-   build_net(df_stage_II_filtered, cor_func = "spearman", n_threads =1)
net_stage_III <-   build_net(df_stage_III_filtered, cor_func = "spearman", n_threads =1)

sub_clusters_modules_stage_I<- get_sub_clusters(net_stage_I$network)
sub_clusters_modules_stage_II<- get_sub_clusters(net_stage_II$network)
sub_clusters_modules_stage_III<- get_sub_clusters(net_stage_III$network)
############################################################################################################################################
# for each module of stage I
for (module in 1:length(unique(sub_clusters_modules_stage_I$sub_module)))
{
  # Copy net and subclusters
  net_stage<-net_stage_I
  sub_clusters_modules_stage<-sub_clusters_modules_stage_I
  
  # Take genes of the module
  module_genes<-sub_clusters_modules_stage[which(sub_clusters_modules_stage$sub_module==module),"gene"]

  if(length(module_genes)>1)
  {
    # build_graph_from_sq_mat use modules genes to crete the graph
    graph_subcluster <- build_graph_from_sq_mat(net_stage$network[module_genes,module_genes])
  
    #######################################################################################################################################
     upper_weight_th = 0.990
    
    # Take subgraph
    cor_subgraph<-net_stage$network[module_genes,module_genes]
    
    # Subclusters
    cor_subgraph[lower.tri(cor_subgraph)] <- NA
    
    # Melt data
    net_stage_subcluster_correlation_network<-melt(cor_subgraph)
    
    net_stage_subcluster_correlation_network<-na.omit(net_stage_subcluster_correlation_network[net_stage_subcluster_correlation_network$value>=upper_weight_th,])
    #######################################################################################################################################
    interactions_stage_subcluster<-unique(net_stage_subcluster_correlation_network[,c(1,2)])
    #######################################################################################################################################
    # If at least one of the genes in the pair are in the interactome
    interactome_data_stage_subcluster<-interactions_stage_subcluster
    
    # Rename columns
    colnames(interactome_data_stage_subcluster)<-c("Gene1","Gene2")
    ########################################################################################################################################
    interactome_data_stage_subcluster<-unique(interactome_data_stage_subcluster[,c("Gene1","Gene2")])
    ########################################################################################################################################
    df_stageI_connectivity   <-unique(data.frame(Conectivity=table(c(interactome_data_stage_subcluster$Gene1,interactome_data_stage_subcluster$Gene2))))
    ########################################################################################################################################
    colnames(df_stageI_connectivity)<-c("Gene","Conectivity")
    ########################################################################################################################################
    # Table for the calculation of entropy
    df_entropy_calulation   <-data.frame(table(df_stageI_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
    
    # Rename colnames
    colnames(df_entropy_calulation)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
    
    # Calculate p(k)
    df_entropy_calulation$p_k<-df_entropy_calulation$count/sum(df_entropy_calulation$count)
    
    # Calculate log2(p(k))
    df_entropy_calulation$log2_pk<-log(df_entropy_calulation$p_k,2)
    
    # Calculate p(k)*log2(p(k))
    df_entropy_calulation$p_k_mult_log2_pk<-df_entropy_calulation$p_k*df_entropy_calulation$log2_pk
    
    # Caclulate entropy value
    Entropy_stage_subcluster_value_Carels  <-abs(sum(df_entropy_calulation$p_k_mult_log2_pk))
    
    # FindClusters_resolution
    png(filename=paste(output_folder,"Network_","graph_Stage_I_Module",module,".png",sep=""), width = 20, height = 20, res=600, units = "cm")
      layout_stage_I <- plot_module(graph_subcluster, upper_weight_th = upper_weight_th,vertex.label.cex = 0.7, node_scaling_max = 7,  legend_cex = 1,title = paste("Network for genes of Stage I\n",paste("Module",module,sep=" "),"\nEntropy : ", round(Entropy_stage_I_value_Carels,4),sep=""))
    dev.off()
  }
}




