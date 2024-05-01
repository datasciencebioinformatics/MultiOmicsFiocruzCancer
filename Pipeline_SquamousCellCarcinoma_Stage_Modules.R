library(igraph)
library(ggnet)
library(network)
library(ggplot2)
###################################################################################################################################################
# Load conversion table
Table1_data<-read.table(file = "/home/felipe/Documentos/LungPortal/EnsemblToUniprotKBconversionList.txt", sep = '\t', header = TRUE,fill=TRUE)    
colnames(Table1_data)<-c("gene_id","gene_symbol")
Table1_data <- Table1_data[match(unique(Table1_data$gene_id), Table1_data$gene_id),]
rownames(Table1_data)<-Table1_data$gene_id
###################################################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output/"   
###################################################################################################################################################
# File path to gene stages
file_genes_Stage_I   <-   paste(output_dir,"df_stageI_interactome",".tsv",sep="")
file_genes_Stage_II   <-  paste(output_dir,"df_stageII_interactome",".tsv",sep="")
file_genes_Stage_III   <- paste(output_dir,"df_stageIII_interactome",".tsv",sep="")

# Gene table
df_stageI_interactome       <-read.table(file = file_genes_Stage_I, sep = '\t', header = TRUE,fill=TRUE)         
df_stageII_interactome      <-read.table(file = file_genes_Stage_II, sep = '\t', header = TRUE,fill=TRUE)
df_stageIII_interactome     <-read.table(file = file_genes_Stage_III, sep = '\t', header = TRUE,fill=TRUE) 

df_stageI_interactome       <-na.omit(data.frame(Gene1=Table1_data[df_stageI_interactome$Gene1,"gene_symbol"],Gene2=Table1_data[df_stageI_interactome$Gene2,"gene_symbol"]))
df_stageII_interactome      <-na.omit(data.frame(Gene1=Table1_data[df_stageII_interactome$Gene1,"gene_symbol"],Gene2=Table1_data[df_stageII_interactome$Gene2,"gene_symbol"]))
df_stageIII_interactome     <-na.omit(data.frame(Gene1=Table1_data[df_stageIII_interactome$Gene1,"gene_symbol"],Gene2=Table1_data[df_stageIII_interactome$Gene2,"gene_symbol"]))

# File path to gene stages
g_stage_I<-graph_from_data_frame(df_stageI_interactome, directed = FALSE, vertices = NULL)
g_stage_II<-graph_from_data_frame(df_stageII_interactome, directed = FALSE, vertices = NULL)
g_stage_III<-graph_from_data_frame(df_stageIII_interactome, directed = FALSE, vertices = NULL)

# Create communities
g_stage_I_ceb <- cluster_optimal(g_stage_I) 
g_stage_I_ceb <- cluster_leiden(g_stage_I_ceb, objective_function = "modularity",n_iterations = 3, resolution_parameter = 0.75)

# Plot communities
df_stageI_interactome.net <- network(df_stageI_interactome, directed = FALSE)
df_stageI_interactome.net %v% "membership" = g_stage_I_ceb$membership
plot_stage_I_labels<-ggnet2(df_stageI_interactome.net,label = g_stage_I_ceb$names, palette = "Set2", color = "membership",layout.par = list(niter = 100000,cell.jitter = 0.75), alpha = 1, size = 10, edge.alpha = 0.2) +  coord_equal()
plot_stage_I_nolabels<-ggnet2(df_stageI_interactome.net, palette = "Set2", color = "membership",layout.par = list(niter = 100000,cell.jitter = 0.75), alpha = 1, size = 10, edge.alpha = 0.2) +  coord_equal()

