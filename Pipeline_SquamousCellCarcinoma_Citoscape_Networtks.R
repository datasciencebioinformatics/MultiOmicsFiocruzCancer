############################################################################################################################################
library(RCy3)
library(RColorBrewer)
# cyrestGET(operation = NULL, parameters = NULL, base.url = "http://localhost:1234")
############################################################################################################################################
stage_I_interaction_file               <- "/media/felipe/UBUNTU\ 20_0/Network_cytoscape/df_stageI_interactome_interactome.tsv" #
stage_II_interaction_file               <- "/media/felipe/UBUNTU\ 20_0/Network_cytoscape/df_stageII_interactome_interactome.tsv" #
stage_III_interaction_file               <- "/media/felipe/UBUNTU\ 20_0/Network_cytoscape/df_stageIII_interactome_interactome.tsv" #
#######################################################################################################################################
stage_I_interaction_data     <-read.table(file = stage_I_interaction_file, sep = '\t', header = TRUE,fill=TRUE)    #
stage_II_interaction_data     <-read.table(file = stage_II_interaction_file, sep = '\t', header = TRUE,fill=TRUE)    #
stage_III_interaction_data     <-read.table(file = stage_III_interaction_file, sep = '\t', header = TRUE,fill=TRUE)    #
#######################################################################################################################################
# Merge data.frame
merge_interactome_gene_symbol <- merge_interactome_gene_symbol[match(unique(merge_interactome_gene_symbol$gene_id), merge_interactome_gene_symbol$gene_id),]

# Assert rownames
rownames(merge_interactome_gene_symbol)<-merge_interactome_gene_symbol$gene_id

# Converted gene_ids
stage_I_interaction_data    <-na.omit(data.frame(Gene1=merge_interactome_gene_symbol[stage_I_interaction_data$Gene1,"gene_symbol"],Gene2=merge_interactome_gene_symbol[stage_I_interaction_data$Gene2,"gene_symbol"]))
stage_II_interaction_data   <-na.omit(data.frame(Gene1=merge_interactome_gene_symbol[stage_II_interaction_data$Gene1,"gene_symbol"],Gene2=merge_interactome_gene_symbol[stage_II_interaction_data$Gene2,"gene_symbol"]))
stage_III_interaction_data  <-na.omit(data.frame(Gene1=merge_interactome_gene_symbol[stage_III_interaction_data$Gene1,"gene_symbol"],Gene2=merge_interactome_gene_symbol[stage_III_interaction_data$Gene2,"gene_symbol"]))

# Read igraph
interactome_igraph_stage_I<-graph_from_data_frame(stage_I_interaction_data, directed = FALSE, vertices = NULL)
interactome_igraph_stage_II<-graph_from_data_frame(stage_II_interaction_data, directed = FALSE, vertices = NULL)
interactome_igraph_stage_III<-graph_from_data_frame(stage_III_interaction_data, directed = FALSE, vertices = NULL)

# greedy method (hiearchical, fast method)
cluster_igraph_stage_I  <- cluster_fast_greedy(interactome_igraph_stage_I)
cluster_igraph_stage_II <- cluster_fast_greedy(interactome_igraph_stage_II)
cluster_igraph_stage_III <- cluster_fast_greedy(interactome_igraph_stage_III)

# Set colour
nb.cols  <- max(c(length(unique(membership(cluster_igraph_stage_I))),length(unique(membership(cluster_igraph_stage_II))), length(unique(membership(cluster_igraph_stage_III)))))
mycolors <- data.frame(colour=colorRampPalette(brewer.pal(8, "Set3"))(nb.cols))

# Generate colors based on media type:
V(interactome_igraph_stage_I)$color  <- mycolors[as.numeric( membership(cluster_igraph_stage_I)),]
V(interactome_igraph_stage_II)$color <- mycolors[as.numeric( membership(cluster_igraph_stage_II)),]
V(interactome_igraph_stage_III)$color<- mycolors[as.numeric( membership(cluster_igraph_stage_III)),]

# Set node size based on audience size:
V(interactome_igraph_stage_I)$size <- degree(interactome_igraph_stage_I)
V(interactome_igraph_stage_II)$size <- degree(interactome_igraph_stage_II)
V(interactome_igraph_stage_III)$size <- degree(interactome_igraph_stage_III)


# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"Panel_subgraph_igraph_stage_I.png",sep=""), width = 30, height = 30, res=1200, units = "cm")                                                                                                    #
  plot(interactome_igraph_stage_I,layout=layout_with_dh, edge.color	="black")
dev.off() 

# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"Panel_subgraph_igraph_stage_II.png",sep=""), width = 30, height = 30, res=1200, units = "cm")                                                                                                    #
  plot(interactome_igraph_stage_II,layout= layout_with_dh, edge.color="black")
dev.off() 

# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"Panel_subgraph_igraph_stage_III.png",sep=""), width = 30, height = 30, res=1200, units = "cm")                                                                                                    #
  plot(interactome_igraph_stage_III,layout=   layout_with_dh, edge.color="black")
dev.off() 


# Save TSV file with genes from Stage1
write_tsv(stage_I_interaction_data, paste(output_dir,"interactome_smbols_stage_I",".tsv",sep=""))
write_tsv(stage_II_interaction_data, paste(output_dir,"interactome_smbols_stage_II",".tsv",sep=""))
write_tsv(stage_III_interaction_data, paste(output_dir,"interactome_smbols_stage_III",".tsv",sep=""))

