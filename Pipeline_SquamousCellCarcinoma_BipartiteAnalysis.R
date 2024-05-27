############################################################################################################################################
library(RCy3)
library(RColorBrewer)
# cyrestGET(operation = NULL, parameters = NULL, base.url = "http://localhost:1234")
# Analyses with the combination of parameter line 119 of the document Parametrization.xlsx
# ≥3	≥1	≤0.05	≥0.85	5456	1798/25	1887/70	1991/182	204/191/1.3396	225/207/1.4054	242/206/1.2978	1276/3819/3.7205	1345/4143/3.7816	1440/4646/3.8299
# /home/felipe/Documentos/scripts_Table7/Script15.R
############################################################################################################################################
interactome_network_Stage_I   <-"/home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/df_stageI_interactome_interactome.tsv"
interactome_network_Stage_II  <-"/home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/df_stageII_interactome_interactome.tsv"
interactome_network_Stage_III <-"/home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/df_stageIII_interactome_interactome.tsv"
##################################################################################################
interactome_network_Stage_I        <-read.table(file = interactome_network_Stage_I, sep = '\t', header = TRUE,fill=TRUE)    #
interactome_network_Stage_II       <-read.table(file = interactome_network_Stage_II, sep = '\t', header = TRUE,fill=TRUE)    #
interactome_network_Stage_III      <-read.table(file = interactome_network_Stage_III, sep = '\t', header = TRUE,fill=TRUE)    #
#######################################################################################################################################
file_unique_gene_stages_I    <-"/home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/DE_GenesPerStageMeansFromPairedUp_unique_stage_I.tsv"
file_unique_gene_stages_II   <-"/home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/DE_GenesPerStageMeansFromPairedUp_unique_stage_II.tsv"
file_unique_gene_stages_III  <-"/home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/DE_GenesPerStageMeansFromPairedUp_unique_stage_III.tsv"

file_gene_stages_I    <-"/home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/DE_GenesPerStageMeansFromPairedUp_Stage_sample_stage_I.tsv"
file_gene_stages_II   <-"/home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/DE_GenesPerStageMeansFromPairedUp_Stage_sample_stage_II.tsv"
file_gene_stages_III  <-"/home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/DE_GenesPerStageMeansFromPairedUp_Stage_sample_stage_III.tsv"

genes_Stage_I                <-read.table(file = file_gene_stages_I, sep = '\t', header = TRUE,fill=TRUE)    #
genes_Stage_II               <-read.table(file = file_gene_stages_II, sep = '\t', header = TRUE,fill=TRUE)    #
genes_Stage_III              <-read.table(file = file_gene_stages_III, sep = '\t', header = TRUE,fill=TRUE)    #

genes_unique_Stage_I         <-read.table(file = file_unique_gene_stages_I, sep = '\t', header = TRUE,fill=TRUE)    #
genes_unique_Stage_II        <-read.table(file = file_unique_gene_stages_II, sep = '\t', header = TRUE,fill=TRUE)    #
genes_unique_Stage_III       <-read.table(file = file_unique_gene_stages_III, sep = '\t', header = TRUE,fill=TRUE)    #
#######################################################################################################################################
# Merge data.frame
merge_interactome_gene_symbol <- merge_interactome_gene_symbol[match(unique(merge_interactome_gene_symbol$gene_id), merge_interactome_gene_symbol$gene_id),]

# Assert rownames
rownames(merge_interactome_gene_symbol)<-merge_interactome_gene_symbol$gene_id

# Converted gene_ids
intractome_network_Stage_I    <-data.frame(Gene1=merge_interactome_gene_symbol[interactome_network_Stage_I$Gene1,"gene_symbol"],Gene2=merge_interactome_gene_symbol[interactome_network_Stage_I$Gene2,"gene_symbol"])
intractome_network_Stage_II   <-data.frame(Gene1=merge_interactome_gene_symbol[interactome_network_Stage_II$Gene1,"gene_symbol"],Gene2=merge_interactome_gene_symbol[interactome_network_Stage_II$Gene2,"gene_symbol"])
intractome_network_Stage_III  <-data.frame(Gene1=merge_interactome_gene_symbol[interactome_network_Stage_III$Gene1,"gene_symbol"],Gene2=merge_interactome_gene_symbol[interactome_network_Stage_III$Gene2,"gene_symbol"])

# Converted gene_ids
symbols_unique_Stage_I        <-merge_interactome_gene_symbol[genes_unique_Stage_I$gene,"gene_symbol"]
symbols_unique_Stage_II       <-merge_interactome_gene_symbol[genes_unique_Stage_II$gene,"gene_symbol"]
symbols_unique_Stage_III      <-merge_interactome_gene_symbol[genes_unique_Stage_III$gene,"gene_symbol"]

# Converted gene_ids
symbols_Stage_I        <-merge_interactome_gene_symbol[genes_Stage_I$gene,"gene_symbol"]
symbols_Stage_II       <-merge_interactome_gene_symbol[genes_Stage_II$gene,"gene_symbol"]
symbols_Stage_III      <-merge_interactome_gene_symbol[genes_Stage_III$gene,"gene_symbol"]
#######################################################################################################################################
# Read igraph
interactome_network_Stage_I<-graph_from_data_frame(interactome_network_Stage_I, directed = FALSE, vertices = NULL)
interactome_network_Stage_II<-graph_from_data_frame(interactome_network_Stage_II, directed = FALSE, vertices = NULL)
interactome_network_Stage_III<-graph_from_data_frame(interactome_network_Stage_III, directed = FALSE, vertices = NULL)

# greedy method (hiearchical, fast method)
cluster_Stage_I = cluster_fast_greedy(interactome_network_Stage_I)
cluster_Stage_II = cluster_fast_greedy(interactome_network_Stage_II)
cluster_Stage_III = cluster_fast_greedy(interactome_network_Stage_III)

########################################################################################################################################
# Calculate gse all stages
# Gene Set Enrichment
# First, I will create go terms per stage using gseGO. 
# Second, I will create dotplots with functional annotation per stage.
# Third, enrichment map will be created.
# Fourth, category netplots and ridgeplot will be created.

# Techinical comments.
# GSE says keytype is not supported. I have tried :
# ENSMBL has been tried and the message remains.
# The identifier with and withour the variant has been tried.
# I am reluctantly to convert the identifiers without testing a valid case. And example:
# I have converted the ids to 101927581,6046,1588,142,3010,57054

# For each gene in stage I
genes_ids_stage_I<-c()
genes_ids_stage_II<-c()
genes_ids_stage_III<-c()
genes_ids_all<-c()

for (gene_id in genes_Stage_I$gene)
{
  # Store gene id in the vector
  genes_ids_stage_I<-c(genes_ids_stage_I,strsplit(gene_id, split = "\\.")[[1]][1])
}
# For each gene in stage II
for (gene_id in genes_Stage_II$gene)
{
  # Store gene id in the vector
  genes_ids_stage_II<-c(genes_ids_stage_II,strsplit(gene_id, split = "\\.")[[1]][1])
}
# For each gene in stage III
for (gene_id in genes_Stage_III$gene)
{
  # Store gene id in the vector
  genes_ids_stage_III<-c(genes_ids_stage_III,strsplit(gene_id, split = "\\.")[[1]][1])
}
# For each gene in stage III
for (gene_id in rownames(unstranded_data_filter))
{
  # Store gene id in the vector
  genes_ids_all<-c(genes_ids_all,strsplit(gene_id, split = "\\.")[[1]][1])
}
########################################################################################################################################
# ids_stage_I
ids_stage_I   <-bitr(genes_ids_stage_I, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
ids_stage_II  <-bitr(genes_ids_stage_II, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
ids_stage_III  <-bitr(genes_ids_stage_III, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
genes_ids_all  <-bitr(genes_ids_all, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")

gse_ALL_Stage_I  <- enrichGO(gene = ids_stage_I$ENTREZID, universe = genes_ids_all,  OrgDb  = org.Hs.eg.db,   ont = "ALL",  pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable = TRUE)
gse_ALL_Stage_II <- enrichGO(gene = ids_stage_II$ENTREZID, universe = genes_ids_all,  OrgDb  = org.Hs.eg.db,   ont = "ALL",  pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable = TRUE)
gse_ALL_Stage_III <- enrichGO(gene = ids_stage_III$ENTREZID, universe = genes_ids_all,  OrgDb  = org.Hs.eg.db,   ont = "ALL",  pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable = TRUE)

gse_MF_Stage_I   <-gse_ALL_Stage_I[gse_ALL_Stage_I$ONTOLOGY=="MF",]
gse_MF_Stage_II  <-gse_ALL_Stage_II[gse_ALL_Stage_II$ONTOLOGY=="MF",]
gse_MF_Stage_III  <-gse_ALL_Stage_III[gse_ALL_Stage_III$ONTOLOGY=="MF",]

gse_BP_Stage_I   <-gse_ALL_Stage_I[gse_ALL_Stage_I$ONTOLOGY=="BP",]
gse_BP_Stage_II  <-gse_ALL_Stage_II[gse_ALL_Stage_II$ONTOLOGY=="BP",]
gse_BP_Stage_III  <-gse_ALL_Stage_III[gse_ALL_Stage_III$ONTOLOGY=="BP",]

gse_CC_Stage_I   <-gse_ALL_Stage_I[gse_ALL_Stage_I$ONTOLOGY=="CC",]
gse_CC_Stage_II  <-gse_ALL_Stage_II[gse_ALL_Stage_II$ONTOLOGY=="CC",]
gse_CC_Stage_III  <-gse_ALL_Stage_III[gse_ALL_Stage_III$ONTOLOGY=="CC",]

gse_MF_Stage_I$Stage<-"Stage I"
gse_MF_Stage_II$Stage<-"Stage II"
gse_MF_Stage_III$Stage<-"Stage III"

gse_BP_Stage_I$Stage<-"Stage I"
gse_BP_Stage_II$Stage<-"Stage II"
gse_BP_Stage_III$Stage<-"Stage III"

gse_CC_Stage_I$Stage<-"Stage I"
gse_CC_Stage_II$Stage<-"Stage II"
gse_CC_Stage_III$Stage<-"Stage III"

# Biological process
gse_BP_Stages<-data.frame(rbind(gse_BP_Stage_I[order(gse_BP_Stage_I$p.adjust),][1:50,],
gse_BP_Stage_II[order(gse_BP_Stage_II$p.adjust),][1:50,],
gse_BP_Stage_III[order(gse_BP_Stage_III$p.adjust),][1:50,]))

# Molecular function
gse_MF_Stages<-data.frame(rbind(gse_MF_Stage_I[order(gse_MF_Stage_I$p.adjust),][1:50,],
gse_MF_Stage_II[order(gse_MF_Stage_II$p.adjust),][1:50,],
gse_MF_Stage_III[order(gse_MF_Stage_III$p.adjust),][1:50,]))

# Celular function  
gse_CC_Stages<-data.frame(rbind(gse_CC_Stage_I[order(gse_CC_Stage_I$p.adjust),][1:50,],
gse_CC_Stage_II[order(gse_CC_Stage_II$p.adjust),][1:50,],
gse_CC_Stage_III[order(gse_CC_Stage_III$p.adjust),][1:50,]))
########################################################################################################################################
plot_bp<-ggplot(gse_BP_Stages, aes(x=Description, y=Count, label=Count)) +geom_bar(stat='identity', aes(fill=p.adjust), width=.5) + coord_flip() + facet_grid(cols = vars(Stage))+ theme_bw() + ggtitle("Biological process")
plot_mf<-ggplot(gse_MF_Stages, aes(x=Description, y=Count, label=Count)) +geom_bar(stat='identity', aes(fill=p.adjust), width=.5) + coord_flip() + facet_grid(cols = vars(Stage))+ theme_bw() + ggtitle("Molecular function")
plot_cc<-ggplot(gse_CC_Stages, aes(x=Description, y=Count, label=Count)) +geom_bar(stat='identity', aes(fill=p.adjust), width=.5) + coord_flip() + facet_grid(cols = vars(Stage))+ theme_bw() + ggtitle("Celular function")
#######################################################################################################################################
# FindClusters_resolution
png(filename=paste(output_folder,"Plot_biological_process.png",sep=""), width = 23, height = 16, res=600, units = "cm")
	plot_bp
dev.off()
# FindClusters_resolution
png(filename=paste(output_folder,"Plot_molecular_function.png",sep=""), width = 23, height = 16, res=600, units = "cm")
	plot_mf
dev.off()
# FindClusters_resolution
png(filename=paste(output_folder,"Plot_ceclular_function.png",sep=""), width = 23, height = 16, res=600, units = "cm")
	plot_cc
dev.off()
#######################################################################################################################################
## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)


t1 <- try(system(paste("rm -r /home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/clusters/",sep=""), intern = TRUE))
t1 <- try(system(paste("mkdir /home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/clusters/",sep=""), intern = TRUE))
########################################################################################################################################
# List of the in stage I
genes_stage_I_per_cluster<-list()

# for each cluster, sabe in file
for (cluster in membership(cluster_Stage_I))
{
  write_tsv(data.frame(Genes=names(which(membership(cluster_Stage_I)==cluster))), paste(output_dir,"/clusters/stage_I_cluster_",cluster,".tsv",sep=""))  
  genes_stage_I_per_cluster[[cluster]]<-names(which(membership(cluster_Stage_I)==cluster))
}
genes_stage_I_annotation <- gseGO(geneList=genes_stage_I_per_cluster,  ont ="ALL",  keyType = "ENSEMBL",  nPerm = 10000,  minGSSize = 3,  maxGSSize = 800,pvalueCutoff = 0.05, verbose = TRUE, OrgDb = organism, pAdjustMethod = "none")
########################################################################################################################################
# List of the in stage II
genes_stage_II_per_cluster<-list()

# for each cluster, sabe in file
for (cluster in membership(cluster_Stage_II))
{
  write_tsv(data.frame(Genes=names(which(membership(cluster_Stage_I)==cluster))), paste(output_dir,"/clusters/stage_II_cluster_",cluster,".tsv",sep=""))
}
########################################################################################################################################
# List of the in stage III
genes_stage_II_per_cluster<-list()

# for each cluster, sabe in file
for (cluster in membership(cluster_Stage_II))
{
  write_tsv(data.frame(Genes=names(which(membership(cluster_Stage_I)==cluster))), paste(output_dir,"/clusters/stage_III_cluster_",cluster,".tsv",sep=""))
}
########################################################################################################################################


########################################################################################################################################
# Colour pallets - R
library("viridis")
viridis(10)

df_colours_stage_I<-data.frame(Colour=viridis(length(unique(membership(cluster_Stage_I)))),Cluster=unique(membership(cluster_Stage_I)))
df_colours_stage_II<-data.frame(Colour=viridis(length(unique(membership(cluster_Stage_II)))),Cluster=unique(membership(cluster_Stage_II)))
df_colours_stage_III<-data.frame(Colour=viridis(length(unique(membership(cluster_Stage_III)))),Cluster=unique(membership(cluster_Stage_III)))

# Set colours according to membership
df_colours_stage_I<-df_colours_stage_I[order(df_colours_stage_I$Cluster),]
df_colours_stage_II<-df_colours_stage_II[order(df_colours_stage_II$Cluster),]
df_colours_stage_III<-df_colours_stage_III[order(df_colours_stage_III$Cluster),]

# Set colours according to membership
rownames(df_colours_stage_I)<-df_colours_stage_I$Cluster
rownames(df_colours_stage_II)<-df_colours_stage_II$Cluster
rownames(df_colours_stage_III)<-df_colours_stage_III$Cluster

V(interactome_network_Stage_I)$color <- df_colours_stage_I[membership(cluster_Stage_I),"Colour"]
V(interactome_network_Stage_II)$color <- df_colours_stage_I[membership(cluster_Stage_II),"Colour"]
V(interactome_network_Stage_III)$color <- df_colours_stage_I[membership(cluster_Stage_III),"Colour"]

# Add n genes per clucster
df_genes_per_stage_I<-data.frame(table(membership(cluster_Stage_I)))
df_genes_per_stage_II<-data.frame(table(membership(cluster_Stage_II)))
df_genes_per_stage_III<-data.frame(table(membership(cluster_Stage_III)))

rownames(df_genes_per_stage_I)<-df_genes_per_stage_I$Var1
rownames(df_genes_per_stage_II)<-df_genes_per_stage_II$Var1
rownames(df_genes_per_stage_III)<-df_genes_per_stage_III$Var1

df_colours_stage_I$Genes<-df_genes_per_stage_I[df_colours_stage_I$Cluster,"Freq"]
df_colours_stage_II$Genes<-df_genes_per_stage_II[df_colours_stage_II$Cluster,"Freq"]
df_colours_stage_III$Genes<-df_genes_per_stage_III[df_colours_stage_III$Cluster,"Freq"]
########################################################################################################################################
# Set number of edges
df_colours_stage_I$nEdges<-0
df_colours_stage_II$nEdges<-0
df_colours_stage_III$nEdges<-0

# Calculate number of edges per cluster
for (cluster_stage_I in rownames(df_colours_stage_I))
{
  # Set number of edges
  df_colours_stage_I[cluster_stage_I,"nEdges"]<-length(E(subgraph(interactome_network_Stage_I, vids=names(membership(cluster_Stage_I)[membership(cluster_Stage_I)==cluster_stage_I]))))
}
# Calculate number of edges per cluster
for (cluster_stage_II in rownames(df_colours_stage_II))
{
  # Set number of edges
  df_colours_stage_II[cluster_stage_II,"nEdges"]<-length(E(subgraph(interactome_network_Stage_II, vids=names(membership(cluster_Stage_II)[membership(cluster_Stage_II)==cluster_stage_II]))))
}
# Calculate number of edges per cluster
for (cluster_stage_III in rownames(df_colours_stage_III))
{
  # Set number of edges
  df_colours_stage_III[cluster_stage_III,"nEdges"]<-length(E(subgraph(interactome_network_Stage_III, vids=names(membership(cluster_Stage_III)[membership(cluster_Stage_III)==cluster_stage_III]))))
}



########################################################################################################################################
V(interactome_network_Stage_I)$size <- log(degree(interactome_network_Stage_I)+0.0001,2)
V(interactome_network_Stage_II)$size <- log(degree(interactome_network_Stage_II)+0.0001,2)
V(interactome_network_Stage_III)$size <- log(degree(interactome_network_Stage_III)+0.0001,2)
########################################################################################################################################
write_tsv(df_colours_stage_I, paste(output_dir,"df_colours_stage_I.tsv",sep=""))
write_tsv(df_colours_stage_II, paste(output_dir,"df_colours_stage_II.tsv",sep=""))
write_tsv(df_colours_stage_III, paste(output_dir,"df_colours_stage_III.tsv",sep=""))


########################################################################################################################################
# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"Panel_subgraph_interactome_stage_I.png",sep=""), width = 30, height = 30, res=600, units = "cm")                                                                                                    #    
  plot(interactome_network_Stage_I,layout=    layout_with_mds, edge.color	="grey50", vertex.label=NA, main="Stage I")          
  legend("right", legend = paste("Cluster ",df_colours_stage_I$Cluster, " :",df_colours_stage_I$Genes,"/",df_colours_stage_I$nEdges,sep=""), pch=21, col=df_colours_stage_I$Colour, pt.bg=df_colours_stage_I$Colour, pt.cex=1, cex=.8, bty="n", ncol=1, title="Nº of genes and edges per cluster")
dev.off() 
# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"Panel_subgraph_interactome_stage_II.png",sep=""), width = 30, height = 30, res=600, units = "cm")                                                                                                    #    
  plot(interactome_network_Stage_II,layout=    layout_with_mds, edge.color	="grey50", vertex.label=NA, main="Stage II")          
  legend("right", legend = paste("Cluster ",df_colours_stage_II$Cluster, " :",df_colours_stage_II$Genes,"/",df_colours_stage_II$nEdges,sep=""), pch=21, col=df_colours_stage_I$Colour, pt.bg=df_colours_stage_II$Colour, pt.cex=1, cex=.8, bty="n", ncol=1, title="Nº of genes and edges per cluster")
dev.off() 
# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"Panel_subgraph_interactome_stage_III.png",sep=""), width = 30, height = 30, res=600, units = "cm")                                                                                                    #    
  plot(interactome_network_Stage_III,layout=    layout_with_mds, edge.color	="grey50", vertex.label=NA, main="Stage III")          
  legend("right", legend = paste("Cluster ",df_colours_stage_III$Cluster, " :",df_colours_stage_III$Genes,"/",df_colours_stage_III$nEdges,sep=""), pch=21, col=df_colours_stage_I$Colour, pt.bg=df_colours_stage_II$Colour, pt.cex=1, cex=.8, bty="n", ncol=1, title="Nº of genes and edges per cluster")
dev.off() 





















########################################################################################################################################
library(ggraph)
library(tidygraph)

data.frame(interactome_network_Stage_I)

# Create graph of highschool friendships
graph <- as_tbl_graph(highschool) |>  mutate(Popularity = centrality_degree(mode = 'in'))

# plot using ggraph
ggraph(graph, layout = 'kk') + 
    geom_edge_fan(aes(alpha = after_stat(index)), show.legend = FALSE) + 
    geom_node_point(aes(size = Popularity)) + 
    facet_edges(~year) + 
    theme_graph(foreground = 'steelblue', fg_text_colour = 'white')


# Bi-partite Co-expression vs. Sub-interactome network for stage I
# Bi-partite Co-expression vs. Sub-interactome network for stage II
# Bi-partite Co-expression vs. Sub-interactome network for stage III
################################################################################################
# Combine genes_unique_Stage_I with genes_unique_Stage_III
unique(c(genes_id_vector_stage_I,genes_id_vector_stage_III))
################################################################################################
# Planejamento igraph networks
# 1- Primeira análise, rede de coexpressão e rede de subinteractomas para cada estágio, juntamente com redes bipartidas.
# Rede co-expressão Estágio I <-> Bipartite co-expressão com sub-interactoma Estágio I <-> Rede sub-interactoma Estágio I
# Rede co-expressão Estágio II <-> Bipartite co-expressão com sub-interactoma Estágio II <-> Rede sub-interactoma Estágio II
# Rede co-expressão Estágio III <-> Bipartite co-expressão com sub-interactoma Estágio II <-> Rede sub-interactoma Estágio III

# 2- Segunda análise, rede para cada estágio, juntamente com redes bipartidas para cada par de estágios
Rede estágio I- Bipartite estagios I e II - Rede estágio II - bipartites II e III - Rede estágio III
