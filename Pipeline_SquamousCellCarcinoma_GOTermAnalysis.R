############################################################################################################################################
# cyrestGET(operation = NULL, parameters = NULL, base.url = "http://localhost:1234")
# Analyses with the combination of parameter line 119 of the document Parametrization.xlsx
# ≥3	≥1	≤0.05	≥0.85	5456	1798/25	1887/70	1991/182	204/191/1.3396	225/207/1.4054	242/206/1.2978	1276/3819/3.7205	1345/4143/3.7816	1440/4646/3.8299
# ENSEMBL ids were converted to ENTREZ ids. enrichGO on org.Hs.eg.db was used (pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05, minGSSize = 3) to anotate 23, 62 and 169 genes from stages I, II and III, respectivelly. Then cnetplot was used to show asociations of genes to top 10 categories.
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
ids_stage_I    <-bitr(genes_unique_Stage_I$gene_id, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
ids_stage_II   <-bitr(genes_unique_Stage_II$gene_id, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
ids_stage_III  <-bitr(genes_unique_Stage_III$gene_id, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
########################################################################################################################################
colnames(ids_stage_I)   <-c("gene_id","ENTREZID")
colnames(ids_stage_II)  <-c("gene_id","ENTREZID")
colnames(ids_stage_III) <-c("gene_id","ENTREZID")
########################################################################################################################################
genes_Stage_I  <-merge(genes_unique_Stage_I,ids_stage_I,by="gene_id")
genes_Stage_II <-merge(genes_unique_Stage_II,ids_stage_II,by="gene_id")
genes_Stage_III<-merge(genes_unique_Stage_III,ids_stage_III,by="gene_id")
genes_Stage_ALL<-rbind(genes_Stage_I,genes_Stage_II,genes_Stage_III)
########################################################################################################################################
go_ALL_Stage_I    <- enrichGO(gene = ids_stage_I$SYMBOL,  OrgDb  = org.Hs.eg.db,      ont = "ALL", pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,readable = TRUE, minGSSize = 3,keyType = "SYMBOL")
go_ALL_Stage_II   <- enrichGO(gene = ids_stage_II$ENTREZ, OrgDb  = org.Hs.eg.db,      ont = "ALL",  pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,readable = TRUE, minGSSize = 3,keyType = "SYMBOL")
go_ALL_Stage_III  <- enrichGO(gene = ids_stage_III$ENTREZ,OrgDb  = org.Hs.eg.db,      ont = "ALL",   pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,readable = TRUE, minGSSize = 3,keyType = "SYMBOL")
go_ALL_Stage      <- enrichGO(gene = genes_Stage_ALL$ENTREZ,OrgDb  = org.Hs.eg.db,    ont = "ALL",   pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,readable = TRUE, minGSSize = 3,keyType = "SYMBOL")
########################################################################################################################################
kegg_ALL_Stage_I    <- enrichKEGG(gene = ids_stage_I$SYMBOL,  organism     = 'hsa',    pvalueCutoff = 0.15,keyType = "SYMBOL")
kegg_ALL_Stage_II   <- enrichKEGG(gene = ids_stage_II$ENTREZ,  organism     = 'hsa',    pvalueCutoff = 0.15,keyType = "SYMBOL")
kegg_ALL_Stage_III  <- enrichKEGG(gene = ids_stage_III$ENTREZ,  organism     = 'hsa',    pvalueCutoff = 0.15,keyType = "SYMBOL")
kegg_ALL_Stage      <- enrichKEGG(gene = genes_Stage_ALL$ENTREZ,  organism     = 'hsa',    pvalueCutoff = 0.15,keyType = "SYMBOL")
########################################################################################################################################
# convert ids
genes_Stage_ALL[genes_Stage_ALL$ENTREZID %in% kegg_ALL_Stage_I@result$geneID ,]



# FindClusters_resolution
png(filename=paste(output_folder,"Plot_cnept_plot_Stage_I.png",sep=""), width = 20, height = 20, res=600, units = "cm")
	cnetplot(gse_ALL_Stage_I, showCategory = 10, layout = "kk") + ggtitle("Stage I") 
dev.off()
# FindClusters_resolution
png(filename=paste(output_folder,"Plot_cnept_plot_Stage_II.png",sep=""), width = 20, height = 20, res=600, units = "cm")
	cnetplot(gse_ALL_Stage_II, showCategory = 5, layout = "kk") + ggtitle("Stage II") 
dev.off()
# FindClusters_resolution
png(filename=paste(output_folder,"Plot_cnept_plot_Stage_III.png",sep=""), width = 20, height = 20, res=600, units = "cm")
	cnetplot(gse_ALL_Stage_III, showCategory = 5, layout = "kk") + ggtitle("Stage III") 
dev.off()

# FindClusters_resolution
png(filename=paste(output_folder,"Plot_cnept_genes_Stage_ALL.png",sep=""), width = 30, height = 30, res=600, units = "cm")
	cnetplot(genes_Stage_ALL, showCategory = 10, layout = "kk",color_gene = "black") + ggtitle("All stages") 
dev.off()

########################################################################################################################################
write.xlsx(x=genes_Stage_I,file=paste(output_dir,"unique_genes",".xlsx",sep=""), sheet="Stage I")
write.xlsx(x=genes_Stage_II,file=paste(output_dir,"unique_genes",".xlsx",sep=""), sheet="Stage II")
write.xlsx(x=genes_Stage_III,file=paste(output_dir,"unique_genes",".xlsx",sep=""), sheet="Stage III")

