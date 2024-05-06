############################################################################################################################################
library(RCy3)
# cyrestGET(operation = NULL, parameters = NULL, base.url = "http://localhost:1234")
############################################################################################################################################
output_folder<-"/home/felipe/Documentos/LungPortal/output/modules/"
############################################################################################################################################
df_correlation_net_stage_I<-data.frame(na.omit(unstranded_data_filter[genes_Stage_I$gene,]))
df_correlation_net_stage_II<-data.frame(na.omit(unstranded_data_filter[genes_Stage_II$gene,]))
df_correlation_net_stage_III<-data.frame(na.omit(unstranded_data_filter[genes_Stage_III$gene,]))
#######################################################################################################################################
rownames(interactome_data_stage_I)
rownames(interactome_data_stage_II)
rownames(interactome_data_stage_III)

# Merge data.frame
merge_interactome_gene_symbol <- merge_interactome_gene_symbol[match(unique(merge_interactome_gene_symbol$gene_id), merge_interactome_gene_symbol$gene_id),]

# Assert rownames
rownames(merge_interactome_gene_symbol)<-merge_interactome_gene_symbol$gene_id

# Converted gene_ids
interactome_smbols_stage_I    <-na.omit(data.frame(Gene1=merge_interactome_gene_symbol[interactome_data_stage_I$Gene1,"gene_symbol"],Gene2=merge_interactome_gene_symbol[interactome_data_stage_I$Gene2,"gene_symbol"]))
interactome_smbols_stage_II   <-na.omit(data.frame(Gene1=merge_interactome_gene_symbol[interactome_data_stage_II$Gene1,"gene_symbol"],Gene2=merge_interactome_gene_symbol[interactome_data_stage_II$Gene2,"gene_symbol"]))
interactome_smbols_stage_III  <-na.omit(data.frame(Gene1=merge_interactome_gene_symbol[interactome_data_stage_III$Gene1,"gene_symbol"],Gene2=merge_interactome_gene_symbol[interactome_data_stage_III$Gene2,"gene_symbol"]))

# Read igraph
interactome_igraph_stage_I<-graph_from_data_frame(interactome_smbols_stage_I, directed = FALSE, vertices = NULL)
interactome_igraph_stage_II<-graph_from_data_frame(interactome_smbols_stage_II, directed = FALSE, vertices = NULL)
interactome_igraph_stage_III<-graph_from_data_frame(interactome_smbols_stage_III, directed = FALSE, vertices = NULL)

# greedy method (hiearchical, fast method)
cluster_igraph_stage_I  <- cluster_fast_greedy(interactome_igraph_stage_I)
cluster_igraph_stage_II <- cluster_fast_greedy(interactome_igraph_stage_II)
cluster_igraph_stage_III <- cluster_fast_greedy(interactome_igraph_stage_III)

# memberships of nodes
membership(cluster_igraph_stage_I)

most_populous_community_stage_I<-data.frame(table(as.vector(membership(cluster_igraph_stage_I))))[1,"Var1"]
most_populous_community_stage_II<-data.frame(table(as.vector(membership(cluster_igraph_stage_II))))[1,"Var1"]
most_populous_community_stage_III<-data.frame(table(as.vector(membership(cluster_igraph_stage_III))))[1,"Var1"]

# Number of verices in the sub-network
data.frame(table(as.vector(membership(cluster_igraph_stage_I))))[1,]
data.frame(table(as.vector(membership(cluster_igraph_stage_II))))[1,]
data.frame(table(as.vector(membership(cluster_igraph_stage_III))))[1,]

# Number of verices in the sub-network
most_populous_module_stage_I   <-data.frame(table(as.vector(membership(cluster_igraph_stage_I))))[2,"Var1"]
most_populous_module_stage_II  <-data.frame(table(as.vector(membership(cluster_igraph_stage_II))))[2,"Var1"]
most_populous_module_stage_III <-data.frame(table(as.vector(membership(cluster_igraph_stage_III))))[2,"Var1"]

# Select sbgraphs
subgraph_igraph_stage_I<-subgraph(interactome_igraph_stage_I, which(membership(cluster_igraph_stage_I)==most_populous_module_stage_I))
subgraph_igraph_stage_II<-subgraph(interactome_igraph_stage_II, which(membership(cluster_igraph_stage_II)==most_populous_module_stage_II))
subgraph_igraph_stage_III<-subgraph(interactome_igraph_stage_III, which(membership(cluster_igraph_stage_III)==most_populous_module_stage_III))

# Set node size based on audience size:
V(subgraph_igraph_stage_I)$size <- degree(subgraph_igraph_stage_I)*2
V(subgraph_igraph_stage_II)$size <- degree(subgraph_igraph_stage_II)*2
V(subgraph_igraph_stage_III)$size <- degree(subgraph_igraph_stage_III)*2

# The labels are currently node IDs.
# Setting them to NA will render no labels:
V(subgraph_igraph_stage_I)$label.color <- "black"
V(subgraph_igraph_stage_II)$label.color <- "black"
V(subgraph_igraph_stage_III)$label.color <- "black"

#change arrow size and edge color:
E(subgraph_igraph_stage_I)$arrow.size <- .2
E(subgraph_igraph_stage_I)$edge.color <- "gray80"

# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"Panel_subgraph_igraph_stage_I.png",sep=""), width = 20, height = 20, res=1200, units = "cm")                                                                                                    #
  plot(subgraph_igraph_stage_I,layout=layout_with_dh, vertex.color="grey")
dev.off() 

# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"Panel_subgraph_igraph_stage_II.png",sep=""), width = 20, height = 20, res=1200, units = "cm")                                                                                                    #
  plot(subgraph_igraph_stage_II,layout= layout_with_dh, vertex.color="gray50")
dev.off() 

# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"Panel_subgraph_igraph_stage_IIi.png",sep=""), width = 20, height = 20, res=1200, units = "cm")                                                                                                    #
  plot(subgraph_igraph_stage_III,layout=   layout_with_dh, vertex.color="gray50")
dev.off() 


# Save TSV file with genes from Stage1
write_tsv(interactome_smbols_stage_I, paste(output_dir,"interactome_smbols_stage_I",".tsv",sep=""))
write_tsv(interactome_smbols_stage_II, paste(output_dir,"interactome_smbols_stage_II",".tsv",sep=""))
write_tsv(interactome_smbols_stage_III, paste(output_dir,"interactome_smbols_stage_III",".tsv",sep=""))



# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"Panel_igraph_stage_I.png",sep=""), width = 20, height = 20, res=1200, units = "cm")                                                                                                    #
  plot(interactome_igraph_stage_I,layout=layout_with_dh, vertex.color="grey")
dev.off() 

# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"Panel_igraph_stage_II.png",sep=""), width = 20, height = 20, res=1200, units = "cm")                                                                                                    #
  plot(interactome_igraph_stage_II,layout= layout_with_dh, vertex.color="gray50")
dev.off() 

# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"Panel_igraph_stage_IIi.png",sep=""), width = 20, height = 20, res=1200, units = "cm")                                                                                                    #
  plot(interactome_igraph_stage_III,layout=   layout_with_dh, vertex.color="gray50")
dev.off() 
