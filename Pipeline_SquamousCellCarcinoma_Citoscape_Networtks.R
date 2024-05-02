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
# Filter table
df_stage_I_filtered <- filter_low_var(t(df_correlation_net_stage_I), pct = 0.75, type = "mean")
df_stage_II_filtered <- filter_low_var(t(df_correlation_net_stage_II), pct = 0.75, type = "mean")
df_stage_III_filtered <- filter_low_var(t(df_correlation_net_stage_III), pct = 0.75, type = "mean")

# calculate network
net_stage_I <-   build_net(df_stage_I_filtered, cor_func = "spearman", n_threads =1)
net_stage_II <-   build_net(df_stage_II_filtered, cor_func = "spearman", n_threads =1)
net_stage_III <-   build_net(df_stage_III_filtered, cor_func = "spearman", n_threads =1)
#######################################################################################################################################
#module_stage_I  <- detect_modules(df_stage_I_filtered,  net_stage_I$network,  detailled_result = TRUE,  merge_threshold = 0.25)
#module_stage_II <- detect_modules(df_stage_II_filtered,  net_stage_II$network,  detailled_result = TRUE,  merge_threshold = 0.25)
#module_stage_II <- detect_modules(df_stage_III_filtered,  net_stage_III$network,  detailled_result = TRUE,  merge_threshold = 0.25)

sub_clusters_modules_stage_I<- get_sub_clusters(net_stage_I$network)
sub_clusters_modules_stage_II<- get_sub_clusters(net_stage_II$network)
sub_clusters_modules_stage_III<- get_sub_clusters(net_stage_III$network)

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
interactions_stage_I<-unique(net_stage_I_correlation_network)
interactions_stage_II<-unique(net_stage_II_correlation_network)
interactions_stage_III<-unique(net_stage_III_correlation_network)
#######################################################################################################################################
# Rename columns
colnames(interactions_stage_I)<-c("Gene1","Gene2","spearman")
colnames(interactions_stage_II)<-c("Gene1","Gene2","spearman")
colnames(interactions_stage_III)<-c("Gene1","Gene2","spearman")
########################################################################################################################################
interactome_data_stage_I<-unique(interactome_data_stage_I)
interactome_data_stage_II<-unique(interactome_data_stage_II)
interactome_data_stage_III<-unique(interactome_data_stage_III)
#######################################################################################################################################
# Create nodes nodes_stage_I
nodes_stage_I   <- data.frame(id=sub_clusters_modules_stage_I$gene,  group=sub_clusters_modules_stage_I$sub_module, stringsAsFactors=FALSE)
nodes_stage_II  <- data.frame(id=sub_clusters_modules_stage_II$gene,  group=sub_clusters_modules_stage_II$sub_module, stringsAsFactors=FALSE)
nodes_stage_III <- data.frame(id=sub_clusters_modules_stage_III$gene,  group=sub_clusters_modules_stage_III$sub_module, stringsAsFactors=FALSE)

# Create nodes nodes_stage_I
edges_stage_I <- data.frame(source=interactions_stage_I$Gene1,  target=interactions_stage_I$Gene2,  weight=interactions_stage_I$spearman,  stringsAsFactors=FALSE)
edges_stage_II <- data.frame(source=interactions_stage_II$Gene1,  target=interactions_stage_II$Gene2,  weight=interactions_stage_II$spearman,  stringsAsFactors=FALSE)
edges_stage_III <- data.frame(source=interactions_stage_III$Gene1,  target=interactions_stage_III$Gene2,  weight=interactions_stage_III$spearman,  stringsAsFactors=FALSE)

# Create networks
stage_I_png_file_name <- paste(output_folder,"/","stage_I_network.png",sep="")
stage_II_png_file_name <- paste(output_folder,"/","stage_II_network.png",sep="")
stage_III_png_file_name <- paste(output_folder,"/","stage_III_network.png",sep="")

# Remove temporary files
if(file.exists(stage_I_png_file_name)){file.remove(stage_I_png_file_name) } 
if(file.exists(stage_II_png_file_name)){file.remove(stage_II_png_file_name) } 
if(file.exists(stage_III_png_file_name)){file.remove(stage_III_png_file_name) } 

#######################################################################################################################################
# Create networks
network_stage_I <- createNetworkFromDataFrames(nodes_stage_I,edges_stage_I, title="Stage I genes")
network_stage_II <- createNetworkFromDataFrames(nodes_stage_II,edges_stage_II, title="Stage II genes")
network_stage_III <- createNetworkFromDataFrames(nodes_stage_III,edges_stage_III, title="Stage III genes")

exportImage(stage_I_png_file_name, type = "png",network=network_stage_I )
exportImage(stage_II_png_file_name, type = "png",network=network_stage_II )
exportImage(stage_III_png_file_name, type = "png",network=network_stage_III )
#######################################################################################################################################

#######################################################################################################################################
