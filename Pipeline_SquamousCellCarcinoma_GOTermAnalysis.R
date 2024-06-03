############################################################################################################################################
# cyrestGET(operation = NULL, parameters = NULL, base.url = "http://localhost:1234")
# Analyses with the combination of parameter line 119 of the document Parametrization.xlsx
# ≥3	≥1	≤0.05	≥0.85	5456	1798/25	1887/70	1991/182	204/191/1.3396	225/207/1.4054	242/206/1.2978	1276/3819/3.7205	1345/4143/3.7816	1440/4646/3.8299
# /home/felipe/Documentos/scripts_Table7/Script15.R
#######################################################################################################################################
file_unique_gene_stages_I    <- paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_I",".tsv",sep="")
file_unique_gene_stages_II   <- paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_II",".tsv",sep="")
file_unique_gene_stages_III  <- paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_III",".tsv",sep="")
#######################################################################################################################################
genes_unique_Stage_I         <-read.table(file = file_unique_gene_stages_I, sep = '\t', header = TRUE,fill=TRUE)    #
genes_unique_Stage_II        <-read.table(file = file_unique_gene_stages_II, sep = '\t', header = TRUE,fill=TRUE)    #
genes_unique_Stage_III       <-read.table(file = file_unique_gene_stages_III, sep = '\t', header = TRUE,fill=TRUE)    #
#######################################################################################################################################
# gene_id
genes_Stage_I$gene_id                   <-""
genes_Stage_II$gene_id                  <-""
genes_Stage_III$gene_id                 <-""
#######################################################################################################################################
# Calculate gse all stages
# For each gene in stage I
genes_ids_stage_I<-c()
genes_ids_stage_II<-c()
genes_ids_stage_III<-c()
genes_ids_all<-c()

# For each gene, add gene_id
for (gene_row in rownames(genes_unique_Stage_I))
{	
	# Store gene id in the vector
	genes_unique_Stage_I[gene_row,"gene_id"]<-strsplit(genes_unique_Stage_I[gene_row,"gene"], split = "\\.")[[1]][1]
	
}
# For each gene, add gene_id
for (gene_row in rownames(genes_unique_Stage_II))
{	
	# Store gene id in the vector
	genes_unique_Stage_II[gene_row,"gene_id"]<-strsplit(genes_unique_Stage_II[gene_row,"gene"], split = "\\.")[[1]][1]	
}
# For each gene, add gene_id
for (gene_row in rownames(genes_unique_Stage_III))
{	
	# Store gene id in the vector
	genes_unique_Stage_III[gene_row,"gene_id"]<-strsplit(genes_unique_Stage_III[gene_row,"gene"], split = "\\.")[[1]][1]	
}
########################################################################################################################################
# ids_stage_I
ids_stage_I    <-bitr(genes_unique_Stage_I$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
ids_stage_II   <-bitr(genes_unique_Stage_II$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
ids_stage_III  <-bitr(genes_unique_Stage_III$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
########################################################################################################################################
colnames(ids_stage_I)   <-c("gene_id","ENTREZID")
colnames(ids_stage_II)  <-c("gene_id","ENTREZID")
colnames(ids_stage_III) <-c("gene_id","ENTREZID")
########################################################################################################################################
genes_Stage_I  <-merge(genes_unique_Stage_I,ids_stage_I,by="gene_id")
genes_Stage_II <-merge(genes_unique_Stage_II,ids_stage_II,by="gene_id")
genes_Stage_III<-merge(genes_unique_Stage_III,ids_stage_III,by="gene_id")
########################################################################################################################################
gse_ALL_Stage_I    <- enrichGO(gene = ids_stage_I$ENTREZ,  OrgDb  = org.Hs.eg.db,      ont = "ALL", pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,readable = TRUE, minGSSize = 3)
gse_ALL_Stage_II   <- enrichGO(gene = ids_stage_II$ENTREZ, OrgDb  = org.Hs.eg.db,     ont = "ALL",  pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,readable = TRUE, minGSSize = 3)
gse_ALL_Stage_III  <- enrichGO(gene = ids_stage_III$ENTREZ,OrgDb  = org.Hs.eg.db,    ont = "ALL",   pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,readable = TRUE, minGSSize = 3)
########################################################################################################################################
# FindClusters_resolution
png(filename=paste(output_folder,"Plot_cnept_plot_Stage_I.png",sep=""), width = 20, height = 20, res=600, units = "cm")
	cnetplot(gse_ALL_Stage_I, showCategory = 10, layout = "kk") + ggtitle("Stage I") 
dev.off()
# FindClusters_resolution
png(filename=paste(output_folder,"Plot_cnept_plot_Stage_II.png",sep=""), width = 20, height = 20, res=600, units = "cm")
	cnetplot(gse_ALL_Stage_II, showCategory = 10, layout = "kk") + ggtitle("Stage II") 
dev.off()
# FindClusters_resolution
png(filename=paste(output_folder,"Plot_cnept_plot_Stage_III.png",sep=""), width = 20, height = 20, res=600, units = "cm")
	cnetplot(gse_ALL_Stage_III, showCategory = 10, layout = "kk") + ggtitle("Stage III") 
dev.off()
########################################################################################################################################
