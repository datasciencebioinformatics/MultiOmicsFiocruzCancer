############################################################################################################################################
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
# gene_id
genes_Stage_I$gene_id                   <-""
genes_Stage_II$gene_id                  <-""
genes_Stage_III$gene_id                 <-""
log2change_tumor_control$gene_id        <-""

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

# For each gene, add gene_id
for (gene_row in rownames(genes_Stage_I))
{	
	# Store gene id in the vector
	genes_Stage_I[gene_row,"gene_id"]<-strsplit(genes_Stage_I[gene_row,"gene"], split = "\\.")[[1]][1]
	
}
# For each gene, add gene_id
for (gene_row in rownames(genes_Stage_II))
{	
	# Store gene id in the vector
	genes_Stage_II[gene_row,"gene_id"]<-strsplit(genes_Stage_II[gene_row,"gene"], split = "\\.")[[1]][1]
	
}
# For each gene, add gene_id
for (gene_row in rownames(genes_Stage_III))
{	
	# Store gene id in the vector
	genes_Stage_III[gene_row,"gene_id"]<-strsplit(genes_Stage_III[gene_row,"gene"], split = "\\.")[[1]][1]
	
}
# For each gene, add gene_id
for (gene_row in rownames(log2change_tumor_control))
{
	# Store gene id in the vector
	log2change_tumor_control[gene_row,"gene_id"]<-strsplit(log2change_tumor_control[gene_row,"gene"], split = "\\.")[[1]][1]
}
########################################################################################################################################

# Create vector Stage_I
vector_stage_I<-genes_Stage_I$log2change
names(vector_stage_I)<-genes_Stage_I$ENTREZID

# Create vector Stage_II
vector_stage_II<-genes_Stage_II$log2change
names(vector_stage_II)<-genes_Stage_II$ENTREZID

# Create vector Stage_III
vector_stage_III<-genes_Stage_III$log2change
names(vector_stage_III)<-genes_Stage_III$ENTREZID

# gse_ALL_Stage
gse_ALL_Stage_I     <- gseKEGG(geneList=vector_stage_I[order(-vector_stage_I)] , organism= 'hsa', nPerm        = 10000,minGSSize    = 3, maxGSSize    = 800,    pvalueCutoff = 0.05,              pAdjustMethod = "BH",               keyType       = "kegg")
gse_ALL_Stage_II    <- gseKEGG(geneList=vector_stage_II[order(-vector_stage_II)] , organism= 'hsa', nPerm        = 10000,minGSSize    = 3, maxGSSize    = 800,    pvalueCutoff = 0.05,              pAdjustMethod = "BH",               keyType       = "kegg")
gse_ALL_Stage_III   <- gseKEGG(geneList=vector_stage_III[order(-vector_stage_III)] , organism= 'hsa', nPerm        = 10000,minGSSize    = 3, maxGSSize    = 800,    pvalueCutoff = 0.05,              pAdjustMethod = "BH",               keyType       = "kegg")

gse_ALL_Stage_I@result$Stage<-"Stage I"
gse_ALL_Stage_II@result$Stage<-"Stage II"
gse_ALL_Stage_III@result$Stage<-"Stage III"

expr_stage_I    <-na.omit(unstranded_data_filter[genes_unique_Stage_I$gene,])
expr_stage_II   <-na.omit(unstranded_data_filter[genes_unique_Stage_II$gene,])
expr_stage_III  <-na.omit(unstranded_data_filter[genes_unique_Stage_III$gene,])

# Data frame for id conversion
df_expr_stage_I<-data.frame(Genes=rownames(expr_stage_I),ENTREZID="",genes_id="")
df_expr_stage_II<-data.frame(Genes=rownames(expr_stage_II),ENTREZID="",genes_id="")
df_expr_stage_III<-data.frame(Genes=rownames(expr_stage_III),ENTREZID="",genes_id="")


# For each gene, add gene_id
for (gene_row in rownames(df_expr_stage_I))
{	
	# Convert genes_id and ENTREZID
	df_expr_stage_I[gene_row,"genes_id"]<-strsplit(df_expr_stage_I[gene_row,"Genes"], split = "\\.")[[1]][1]
	try(df_expr_stage_I[gene_row,"ENTREZID"]<-bitr(df_expr_stage_I[gene_row,"genes_id"], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")[1,"ENTREZID"], silent = TRUE)
}
df_expr_stage_I<-na.omit(df_expr_stage_I)

expr_stage_I<-expr_stage_I[df_expr_stage_I$Genes,]
rownames(expr_stage_I)<-df_expr_stage_I$ENTREZID


# Run GSEA 
# First, all msigdbr
pathwaysDF <- msigdbr("human")

# Split name of pathways
pathways <- split(as.character(pathwaysDF$entrez_gene), pathwaysDF$gs_name)

# Run fast gsea 
gesecaRes <- geseca(pathways, expr_stage_I, minSize = 5, maxSize = 500)

plotGesecaTable(gesecaRes , pathways, E=expr_stage_I) + theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
  axis.ticks = element_blank())




# Biological process
gse_KEGG<-na.omit(data.frame(rbind(gse_ALL_Stage_I@result,gse_ALL_Stage_II@result, gse_ALL_Stage_III@result)))
########################################################################################################################################
plot_KEGG<-ggplot(gse_KEGG, aes(x=Description, y=setSize, label=setSize)) +geom_bar(stat='identity', aes(fill=p.adjust), width=.5) + coord_flip() + facet_grid(cols = vars(Stage))+ theme_bw() + ggtitle("KEGG pathway")
#######################################################################################################################################
# FindClusters_resolution
png(filename=paste(output_folder,"Plot_KEGG.png",sep=""), width = 23, height = 16, res=600, units = "cm")
	plot_KEGG
dev.off()
#######################################################################################################################################
write.xlsx(na.omit(gse_KEGG), "KEGG_PATHWAYS", file=paste(output_dir,"/clusters/kegg_pathways.xlsx",sep=""),append = FALSE) # where x is a data.frame with a Date column.
#######################################################################################################################################
