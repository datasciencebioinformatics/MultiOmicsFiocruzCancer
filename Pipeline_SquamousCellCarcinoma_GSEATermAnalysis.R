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
expr_stage_I    <-na.omit(unstranded_data_filter[genes_unique_Stage_I$gene,])
expr_stage_II   <-na.omit(unstranded_data_filter[genes_unique_Stage_II$gene,])
expr_stage_III  <-na.omit(unstranded_data_filter[genes_unique_Stage_III$gene,])

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
#######################################################################################################################################
# Here I must check if if I use all the msigdbr databases or any in particulart
# Run GSEA 
# First, all msigdbr
pathwaysDF <- msigdbr("human")

# Split name of pathways
pathways <- split(as.character(pathwaysDF$entrez_gene), pathwaysDF$gs_name)

# Run fast gsea 
gesecaRes <- geseca(pathways, expr_stage_I, minSize = 5, maxSize = 500)
#######################################################################################################################################
# Here I must check 
plotGesecaTable(gesecaRes , pathways, E=expr_stage_I) + theme(axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank())
#######################################################################################################################################
# FindClusters_resolution
png(filename=paste(output_folder,"Plot_KEGG.png",sep=""), width = 23, height = 16, res=600, units = "cm")
	plot_KEGG
dev.off()
#######################################################################################################################################
write.xlsx(na.omit(gse_KEGG), "KEGG_PATHWAYS", file=paste(output_dir,"/clusters/kegg_pathways.xlsx",sep=""),append = FALSE) # where x is a data.frame with a Date column.
#######################################################################################################################################
