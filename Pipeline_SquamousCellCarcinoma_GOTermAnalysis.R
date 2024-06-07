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
colnames(ids_stage_I)   <-c("gene_id","ENTREZID","SYMBOL")
colnames(ids_stage_II)  <-c("gene_id","ENTREZID","SYMBOL")
colnames(ids_stage_III) <-c("gene_id","ENTREZID","SYMBOL")
########################################################################################################################################
genes_Stage_I  <-merge(genes_unique_Stage_I,ids_stage_I,by="gene_id")
genes_Stage_II <-merge(genes_unique_Stage_II,ids_stage_II,by="gene_id")
genes_Stage_III<-merge(genes_unique_Stage_III,ids_stage_III,by="gene_id")
genes_Stage_ALL<-unique(rbind(ids_stage_I,ids_stage_II,ids_stage_III))
rownames(genes_Stage_ALL)<-genes_Stage_ALL$ENTREZID
########################################################################################################################################
go_ALL_Stage_I    <- enrichGO(gene = ids_stage_I$SYMBOL,  OrgDb  = org.Hs.eg.db,      ont = "ALL", pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,readable = TRUE, minGSSize = 3,keyType = "SYMBOL")
go_ALL_Stage_II   <- enrichGO(gene = ids_stage_II$SYMBOL, OrgDb  = org.Hs.eg.db,      ont = "ALL",  pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,readable = TRUE, minGSSize = 3,keyType = "SYMBOL")
go_ALL_Stage_III  <- enrichGO(gene = ids_stage_III$SYMBOL,OrgDb  = org.Hs.eg.db,      ont = "ALL",   pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,readable = TRUE, minGSSize = 3,keyType = "SYMBOL")
go_ALL_Stage      <- enrichGO(gene = genes_Stage_ALL$SYMBOL,OrgDb  = org.Hs.eg.db,    ont = "ALL",   pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,readable = TRUE, minGSSize = 3,keyType = "SYMBOL")
########################################################################################################################################
kegg_ALL_Stage_I    <- enrichKEGG(gene = ids_stage_I$ENTREZ,  organism     = 'hsa',    pvalueCutoff = 0.15)
kegg_ALL_Stage_II   <- enrichKEGG(gene = ids_stage_II$ENTREZ,  organism     = 'hsa',    pvalueCutoff = 0.15)
kegg_ALL_Stage_III  <- enrichKEGG(gene = ids_stage_III$ENTREZ,  organism     = 'hsa',    pvalueCutoff = 0.15)
kegg_ALL_Stage      <- enrichKEGG(gene = genes_Stage_ALL$ENTREZ,  organism     = 'hsa',    pvalueCutoff = 0.05)
########################################################################################################################################
reactome_ALL_Stage_I    <- enrichPathway(gene=ids_stage_I$ENTREZ, pvalueCutoff = 0.20, readable=TRUE)
reactome_ALL_Stage_II   <- enrichPathway(gene=ids_stage_II$ENTREZ, pvalueCutoff = 0.20, readable=TRUE)
reactome_ALL_Stage_III  <- enrichPathway(gene=ids_stage_III$ENTREZ, pvalueCutoff = 0.20, readable=TRUE)
reactome_ALL_Stage      <- enrichPathway(gene = genes_Stage_ALL$ENTREZ,pvalueCutoff = 0.05, readable=TRUE)
########################################################################################################################################
# A table with the associate go term per gene.
# First, take the results per stage
# Go terms
go_results_Stage_I<-go_ALL_Stage_I@result
go_results_Stage_II<-go_ALL_Stage_II@result
go_results_Stage_III<-go_ALL_Stage_III@result

# Kegg terms
kegg_results_Stage_I<-kegg_ALL_Stage_I@result
kegg_results_Stage_II<-kegg_ALL_Stage_II@result
kegg_results_Stage_III<-kegg_ALL_Stage_III@result

# Reactome terms
reactome_results_Stage_I<-reactome_ALL_Stage_I@result
reactome_results_Stage_II<-reactome_ALL_Stage_II@result
reactome_results_Stage_III<-reactome_ALL_Stage_III@result

# Second, add information about stage
go_results_Stage_I$Stage<-"Stage I"
go_results_Stage_II$Stage<-"Stage II"
go_results_Stage_III$Stage<-"Stage III"

# Second, add information about stage
kegg_results_Stage_I$Stage<-"Stage I"
kegg_results_Stage_II$Stage<-"Stage II"
kegg_results_Stage_III$Stage<-"Stage III"

# Second, add information about stage
reactome_results_Stage_I$Stage<-"Stage I"
reactome_results_Stage_II$Stage<-"Stage II"
reactome_results_Stage_III$Stage<-"Stage III"

# Merge all stages
go_results_all_Stages         <-rbind(rbind(go_results_Stage_I,go_results_Stage_II,go_results_Stage_III)[,c("ID","p.adjust","Description","Count","Stage")]
kegg_results_all_Stages       <-rbind(kegg_results_Stage_I,kegg_results_Stage_II,kegg_results_Stage_III)[,c("ID","p.adjust","Description","Count","Stage")]
reactome_results_all_Stages   <-rbind(reactome_results_Stage_I,reactome_results_Stage_II,reactome_results_Stage_III)[,c("ID","p.adjust","Description","Count","Stage")]

# Merge all stages
go_results_all_Stages$Layer<-"GO"
kegg_results_all_Stages$Layer<-"KEGG"
reactome_results_all_Stages$Layer<-"Reactome"

# Merge all tables
all_anotation_results<-rbind(go_results_all_Stages,kegg_results_all_Stages,reactome_results_all_Stages)
######################################################################################################################
for rowrownames(go_results_all_Stages)









# FindClusters_resolution
png(filename=paste(output_folder,"Plot_KEGG_all_Stage.png",sep=""), width = 45, height = 15, res=600, units = "cm")
	p1<-cnetplot(kegg_ALL_Stage_I, showCategory = 3, layout = "kk",color.params = list(gene =c("black"))) + ggtitle("KEGG pathway - Stage I") 
	p2<-cnetplot(kegg_ALL_Stage_II, showCategory = 3, layout = "kk",color.params = list(gene =c("black"))) + ggtitle("KEGG pathway - Stage II") 
	p3<-cnetplot(kegg_ALL_Stage_III, showCategory = 3, layout = "kk",color.params = list(gene =c("black"))) + ggtitle("KEGG pathway - Stage III") 
	grid.arrange(p1, p2,p3, nrow = 1)	
dev.off()

# FindClusters_resolution
png(filename=paste(output_folder,"Plot_GO_all_Stage.png",sep=""), width = 45, height = 15, res=600, units = "cm")
	p1<-cnetplot(go_ALL_Stage_I, showCategory = 3, layout = "kk",color.params = list(gene =c("black"))) + ggtitle("GO terms - Stage I") 
	p2<-cnetplot(go_ALL_Stage_II, showCategory = 3, layout = "kk",color.params = list(gene =c("black"))) + ggtitle("GO terms - Stage II") 
	p3<-cnetplot(go_ALL_Stage_III, showCategory = 3, layout = "kk",color.params = list(gene =c("black"))) + ggtitle("GO terms - Stage III") 
	grid.arrange(p1, p2,p3, nrow = 1)	
dev.off()

# FindClusters_resolution
png(filename=paste(output_folder,"Plot_Reactome_all_Stage.png",sep=""), width = 45, height = 15, res=600, units = "cm")
	p1<-cnetplot(reactome_ALL_Stage_I, showCategory = 3, layout = "kk",color.params = list(gene =c("black"))) + ggtitle("Reactome pathway - Stage I") 
	p2<-cnetplot(reactome_ALL_Stage_II, showCategory = 3, layout = "kk",color.params = list(gene =c("black"))) + ggtitle("Reactome pathway - Stage II") 
	p3<-cnetplot(reactome_ALL_Stage_III, showCategory = 3, layout = "kk",color.params = list(gene =c("black"))) + ggtitle("Reactome pathway - Stage III") 
	grid.arrange(p1, p2,p3, nrow = 1)	
dev.off()

########################################################################################################################################
write.xlsx(x=genes_Stage_I,file=paste(output_dir,"unique_genes",".xlsx",sep=""), sheet="Stage I")
write.xlsx(x=genes_Stage_II,file=paste(output_dir,"unique_genes",".xlsx",sep=""), sheet="Stage II")
write.xlsx(x=genes_Stage_III,file=paste(output_dir,"unique_genes",".xlsx",sep=""), sheet="Stage III")
