############################################################################################################################################
# cyrestGET(operation = NULL, parameters = NULL, base.url = "http://localhost:1234")
# Analyses with the combination of parameter line 119 of the document Parametrization.xlsx
# ≥3	≥1	≤0.05	≥0.85	5456	1798/25	1887/70	1991/182	204/191/1.3396	225/207/1.4054	242/206/1.2978	1276/3819/3.7205	1345/4143/3.7816	1440/4646/3.8299
# /home/felipe/Documentos/scripts_Table7/Script15.R
############################################################################################################################################
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
# Here I must check what set of genes to use in the GSEA
# if the complete set of genes, or if only the selected set of genes
expr_stage_I    <-unstranded_rpkm[genes_unique_Stage_I$gene,]
expr_stage_II   <-unstranded_rpkm[genes_unique_Stage_II$gene,]
expr_stage_III  <-unstranded_rpkm[genes_unique_Stage_III$gene,]

# Data frame for id conversion
df_expr_stage_I<-data.frame(Genes=rownames(expr_stage_I),ENTREZID="",genes_id="")
df_expr_stage_II<-data.frame(Genes=rownames(expr_stage_II),ENTREZID="",genes_id="")
df_expr_stage_III<-data.frame(Genes=rownames(expr_stage_III),ENTREZID="",genes_id="")
#######################################################################################################################################
# Here I must check if the conversion is convering all the genes.
# For each gene, add gene_id
for (gene_row in rownames(df_expr_stage_I))
{	
	# Convert genes_id and ENTREZID
	df_expr_stage_I[gene_row,"genes_id"]<-strsplit(df_expr_stage_I[gene_row,"Genes"], split = "\\.")[[1]][1]
	try(df_expr_stage_I[gene_row,"ENTREZID"]<-bitr(df_expr_stage_I[gene_row,"genes_id"], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")[1,"ENTREZID"], silent = TRUE)
}
# Set rownames
rownames(df_expr_stage_I)<-df_expr_stage_I$Genes

# Keep only first occcurance
df_expr_stage_I <- df_expr_stage_I[match(unique(df_expr_stage_I$genes_id), df_expr_stage_I$genes_id),]

# Filter dataset
expr_stage_I<-expr_stage_I[df_expr_stage_I$Genes,]

# Set rownames on expr_stage_I
rownames(expr_stage_I)<-df_expr_stage_I[rownames(expr_stage_I),"ENTREZID"]
#######################################################################################################################################
# Here I must check if the conversion is convering all the genes.
# For each gene, add gene_id
for (gene_row in rownames(df_expr_stage_II))
{	
	# Convert genes_id and ENTREZID
	df_expr_stage_II[gene_row,"genes_id"]<-strsplit(df_expr_stage_II[gene_row,"Genes"], split = "\\.")[[1]][1]
	try(df_expr_stage_II[gene_row,"ENTREZID"]<-bitr(df_expr_stage_II[gene_row,"genes_id"], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")[1,"ENTREZID"], silent = TRUE)
}
# Set rownames
rownames(df_expr_stage_II)<-df_expr_stage_II$Genes

# Keep only first occcurance
df_expr_stage_II <- df_expr_stage_II[match(unique(df_expr_stage_II$genes_id), df_expr_stage_II$genes_id),]

# Filter dataset
expr_stage_II<-expr_stage_II[df_expr_stage_II$Genes,]

# Set rownames on expr_stage_I
rownames(expr_stage_II)<-df_expr_stage_II[rownames(expr_stage_II),"ENTREZID"]
#######################################################################################################################################
# Here I must check if the conversion is convering all the genes.
# For each gene, add gene_id
for (gene_row in rownames(df_expr_stage_III))
{	
	# Convert genes_id and ENTREZID
	df_expr_stage_III[gene_row,"genes_id"]<-strsplit(df_expr_stage_III[gene_row,"Genes"], split = "\\.")[[1]][1]
	try(df_expr_stage_III[gene_row,"ENTREZID"]<-bitr(df_expr_stage_III[gene_row,"genes_id"], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")[1,"ENTREZID"], silent = TRUE)
}
# Set rownames
rownames(df_expr_stage_III)<-df_expr_stage_III$Genes

# Keep only first occcurance
df_expr_stage_III <- df_expr_stage_III[match(unique(df_expr_stage_III$genes_id), df_expr_stage_III$genes_id),]

# Filter dataset
expr_stage_III<-expr_stage_III[df_expr_stage_III$Genes,]

# Set rownames on expr_stage_I
rownames(expr_stage_III)<-df_expr_stage_III[rownames(expr_stage_III),"ENTREZID"]
#######################################################################################################################################
# Here I must check if if I use all the msigdbr databases or any in particulart
# Run GSEA 
# First, all msigdbr
pathwaysDF <- msigdbr("human")

# Split name of pathways
pathways <- split(as.character(pathwaysDF$entrez_gene), pathwaysDF$gs_name)
#######################################################################################################################################
# Split name of pathways
stage_I_l2fc<-genes_unique_Stage_I$log2change
names(stage_I_l2fc)<-df_expr_stage_I[genes_unique_Stage_I$gene,"ENTREZID"]

# Split name of pathways
stage_II_l2fc<-genes_unique_Stage_II$log2change
names(stage_II_l2fc)<-df_expr_stage_II[genes_unique_Stage_II$gene,"ENTREZID"]

# Split name of pathways
stage_III_l2fc<-genes_unique_Stage_III$log2change
names(stage_III_l2fc)<-df_expr_stage_III[genes_unique_Stage_II$gene,"ENTREZID"]
#######################################################################################################################################

# Calculate fgsea
fgsaRes_stage_I    <- fgsea(pathways, stage_I_l2fc, minSize=5, maxSize=500)

plotGseaTable(pathways, stage_I_l2fc,fgsaRes_stage_I, gseaParam=0.5)


#######################################################################################################################################
# Here I must check 
plotGseaTable(gesecaRes_stage_I, pathways[gesecaRes_stage_I$pathway], E=expr_stage_I) + theme(axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank())
plotGseaTable(gesecaRes_stage_II, pathways[gesecaRes_stage_II$pathway], E=expr_stage_II) + theme(axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank())
plotGseaTable(gesecaRes_stage_III, pathways[gesecaRes_stage_III$pathway], E=expr_stage_III) + theme(axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank())
#######################################################################################################################################
# FindClusters_resolution
png(filename=paste(output_folder,"Plot_KEGG.png",sep=""), width = 23, height = 16, res=600, units = "cm")
	plot_KEGG
dev.off()
#######################################################################################################################################
write.xlsx(na.omit(gse_KEGG), "KEGG_PATHWAYS", file=paste(output_dir,"/clusters/kegg_pathways.xlsx",sep=""),append = FALSE) # where x is a data.frame with a Date column.
#######################################################################################################################################
