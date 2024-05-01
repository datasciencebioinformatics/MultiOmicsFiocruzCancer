library(igraph)
library(ggnet)
library(network)

# File path to gene stages
file_genes_Stage_I   <-   paste(output_dir,"df_stageI_interactome",".tsv",sep="")
file_genes_Stage_II   <-  paste(output_dir,"df_stageII_interactome",".tsv",sep="")
file_genes_Stage_III   <- paste(output_dir,"df_stageIII_interactome",".tsv",sep="")

# Gene table
df_stageI_interactome       <-read.table(file = file_genes_Stage_I, sep = '\t', header = TRUE,fill=TRUE)         
df_stageII_interactome      <-read.table(file = file_genes_Stage_II, sep = '\t', header = TRUE,fill=TRUE)
df_stageIII_interactome     <-read.table(file = file_genes_Stage_III, sep = '\t', header = TRUE,fill=TRUE) 

# File path to gene stages
g_stage_I<-graph_from_data_frame(df_stageI_interactome, directed = FALSE, vertices = NULL)
g_stage_II<-graph_from_data_frame(df_stageII_interactome, directed = FALSE, vertices = NULL)
g_stage_III<-graph_from_data_frame(df_stageIII_interactome, directed = FALSE, vertices = NULL)



g_stage_I_ceb <- cluster_edge_betweenness(g_stage_I) 
c_stage_I$names[c_stage_I$membership==1]

# greedy method (hiearchical, fast method)
c_stage_I = cluster_fast_greedy(g_stage_I)
coords = layout_with_fr(g_stage_I)
plot(c_stage_I, g_stage_I, layout=coords) 


df_stageI_interactome.net <- network(df_stageI_interactome, directed = FALSE)
df_stageI_interactome.net %v% "membership" = g_stage_I_ceb$membership
ggnet2(df_stageI_interactome.net, color = "membership",layout.par = list(cell.jitter = 0.75)) +  coord_equal()

ggnet2(df_stageI_interactome.net, color = membership) +  coord_equal()


net %v% "x" = x[, 1]

for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(g_stage_I)) 
  plot(g_stage_I, edge.arrow.mode=0, layout=l, main=layout) }
