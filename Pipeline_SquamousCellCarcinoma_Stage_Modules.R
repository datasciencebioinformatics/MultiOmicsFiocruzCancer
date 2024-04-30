library(igraph)

# File path to gene stages
file_genes_Stage_I   <-   paste(output_dir,"df_stageI_interactome",".tsv",sep="")
file_genes_Stage_II   <-  paste(output_dir,"df_stageII_interactome",".tsv",sep="")
file_genes_Stage_III   <- paste(output_dir,"df_stageIII_interactome",".tsv",sep="")

# Gene table
df_stageI_interactome       <-read.table(file = file_genes_Stage_I, sep = '\t', header = TRUE,fill=TRUE)         
df_stageII_interactome      <-read.table(file = file_genes_Stage_II, sep = '\t', header = TRUE,fill=TRUE)
df_stageIII_interactome     <-read.table(file = file_genes_Stage_III, sep = '\t', header = TRUE,fill=TRUE) 

# File path to gene stages
g_stage_I<-graph_from_data_frame(df_stageI_interactome, directed = TRUE, vertices = NULL)
g_stage_II<-graph_from_data_frame(df_stageII_interactome, directed = TRUE, vertices = NULL)
g_stage_III<-graph_from_data_frame(df_stageIII_interactome, directed = TRUE, vertices = NULL)

g_stage_I_ceb <- cluster_edge_betweenness(g_stage_I) 
dendPlot(g_stage_I_ceb, mode="hclust")
l <- do.call("layout_with_kk", list(g_stage_I)) 
plot(g_stage_I_ceb, g_stage_I, layout=l) 

# FindClusters_resolution
png(filename=paste(output_dir,"Network_","Stage_I",".png",sep=""), width = 32, height = 32, res=600, units = "cm")
  plot(g_stage_I_ceb, g_stage_I, layout=l) 
dev.off()


layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
par(mfrow=c(3,3), mar=c(1,1,1,1))

for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(g_stage_I)) 
  plot(g_stage_I, edge.arrow.mode=0, layout=l, main=layout) }
