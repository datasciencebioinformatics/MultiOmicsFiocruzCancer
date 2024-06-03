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
#######################################################################################################################################
# Here I must check what set of genes to use in the GSEA
# if the complete set of genes, or if only the selected set of genes
expr_stage_I    <-unstranded_rpkm[genes_unique_Stage_I$gene,sample_stage_I]
expr_stage_II   <-unstranded_rpkm[genes_unique_Stage_II$gene,sample_stage_II]
expr_stage_III  <-unstranded_rpkm[genes_unique_Stage_III$gene,sample_stage_III]
expr_all        <-unstranded_rpkm

# Data frame for id conversion
df_expr_stage_I<-data.frame(Genes=rownames(expr_stage_I),ENTREZID="",genes_id="")
df_expr_stage_II<-data.frame(Genes=rownames(expr_stage_II),ENTREZID="",genes_id="")
df_expr_stage_III<-data.frame(Genes=rownames(expr_stage_III),ENTREZID="",genes_id="")
df_expr_all<-data.frame(Genes=rownames(unstranded_data_filter),ENTREZID="",genes_id="")
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
df_expr_stage_I <- df_expr_stage_I[match(unique(df_expr_stage_I$ENTREZID), df_expr_stage_I$ENTREZID),]

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
df_expr_stage_II <- df_expr_stage_II[match(unique(df_expr_stage_II$ENTREZID), df_expr_stage_II$ENTREZID),]

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
df_expr_stage_III <- df_expr_stage_III[match(unique(df_expr_stage_III$ENTREZID), df_expr_stage_III$ENTREZID),]

# Filter dataset
expr_stage_III<-expr_stage_III[df_expr_stage_III$Genes,]

# Set rownames on expr_stage_I
rownames(expr_stage_III)<-df_expr_stage_III[rownames(expr_stage_III),"ENTREZID"]
#######################################################################################################################################
# Here I must check if the conversion is convering all the genes.
# For each gene, add gene_id
for (gene_row in rownames(unstranded_data_filter))
{	
	# Convert genes_id and ENTREZID
	df_expr_all[gene_row,"genes_id"]<-strsplit(df_expr_all[gene_row,"Genes"], split = "\\.")[[1]][1]
	try(df_expr_all[gene_row,"ENTREZID"]<-bitr(df_expr_all[gene_row,"genes_id"], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")[1,"ENTREZID"], silent = TRUE)
}
# Set rownames
rownames(df_expr_all)<-df_expr_all$Genes

# Keep only first occcurance
df_expr_all <- df_expr_all[match(unique(df_expr_all$genes_id), df_expr_all$genes_id),]
df_expr_all <- df_expr_all[match(unique(df_expr_all$ENTREZID), df_expr_all$ENTREZID),]

# Filter dataset
expr_all<-expr_all[df_expr_all$Genes,]

# Set rownames on expr_stage_I
rownames(expr_all)<-df_expr_all[rownames(expr_all),"ENTREZID"]

#######################################################################################################################################
# Here I must check if if I use all the msigdbr databases or any in particulart
# Run GSEA 
#######################################################################################################################################
# CC
clusterProfiler_CC_Stage_I   <- groupGO(gene     = df_expr_stage_I$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "CC", level    = 2,   readable = TRUE)@result
clusterProfiler_CC_Stage_II  <- groupGO(gene     = df_expr_stage_II$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "CC", level    = 2,   readable = TRUE)@result
clusterProfiler_CC_Stage_III <- groupGO(gene     = df_expr_stage_III$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "CC", level    = 2,   readable = TRUE)@result

# MF
clusterProfiler_MF_Stage_I   <- groupGO(gene     = df_expr_stage_I$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "MF", level    = 2,   readable = TRUE)@result
clusterProfiler_MF_Stage_II  <- groupGO(gene     = df_expr_stage_II$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "MF", level    = 2,   readable = TRUE)@result
clusterProfiler_MF_Stage_III <- groupGO(gene     = df_expr_stage_III$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "MF", level    = 2,   readable = TRUE)@result

# BP
clusterProfiler_BP_Stage_I   <- groupGO(gene     = df_expr_stage_I$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "BP", level    = 2,   readable = TRUE)@result
clusterProfiler_BP_Stage_II  <- groupGO(gene     = df_expr_stage_II$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "BP", level    = 2,   readable = TRUE)@result
clusterProfiler_BP_Stage_III <- groupGO(gene     = df_expr_stage_III$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "BP", level    = 2,   readable = TRUE)@result
#######################################################################################################################################
# Biological process
clusterProfiler_CC_Stage_I$Stage<-"Stage I"
clusterProfiler_CC_Stage_II$Stage<-"Stage II"
clusterProfiler_CC_Stage_III$Stage<-"Stage III"
clusterProfiler_CC_Stage<-na.omit(data.frame(rbind(clusterProfiler_CC_Stage_I,clusterProfiler_CC_Stage_II,clusterProfiler_CC_Stage_III)))

# A plot with celular component
plot_CC<-ggplot(clusterProfiler_CC_Stage, aes(x=Description, y=Count, label=Count)) +geom_bar(stat='identity', width=.5) + coord_flip() + facet_grid(cols = vars(Stage))+ theme_bw() + ggtitle("Celular component")
#######################################################################################################################################
# Biological process
clusterProfiler_BP_Stage_I$Stage<-"Stage I"
clusterProfiler_BP_Stage_II$Stage<-"Stage II"
clusterProfiler_BP_Stage_III$Stage<-"Stage III"
clusterProfiler_BP_Stage<-na.omit(data.frame(rbind(clusterProfiler_BP_Stage_I,clusterProfiler_BP_Stage_II,clusterProfiler_BP_Stage_III)))

# A plot with celular component
plot_BP<-ggplot(clusterProfiler_BP_Stage, aes(x=Description, y=Count, label=Count)) +geom_bar(stat='identity', width=.5) + coord_flip() + facet_grid(cols = vars(Stage))+ theme_bw() + ggtitle("Biological Process")
#######################################################################################################################################
# Biological process
clusterProfiler_MF_Stage_I$Stage<-"Stage I"
clusterProfiler_MF_Stage_II$Stage<-"Stage II"
clusterProfiler_MF_Stage_III$Stage<-"Stage III"
clusterProfiler_MF_Stage<-na.omit(data.frame(rbind(clusterProfiler_MF_Stage_I,clusterProfiler_MF_Stage_II,clusterProfiler_MF_Stage_III)))

# A plot with celular component
plot_MF<-ggplot(clusterProfiler_MF_Stage, aes(x=Description, y=Count, label=Count)) +geom_bar(stat='identity',  width=.5) + coord_flip() + facet_grid(cols = vars(Stage))+ theme_bw() + ggtitle("Molecular function")
#######################################################################################################################################
# GroupGO
# CC
clusterProfiler_CC_Stage_I   <- enrichGO(gene     = df_expr_stage_I$ENTREZID, universe=, OrgDb    = org.Hs.eg.db, pAdjustMethod = "BH",   pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,  readable = TRUE))@result
clusterProfiler_CC_Stage_II  <- enrichGO(gene     = df_expr_stage_II$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "CC", level    = 2,   readable = TRUE)@result
clusterProfiler_CC_Stage_III <- enrichGO(gene     = df_expr_stage_III$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "CC", level    = 2,   readable = TRUE)@result

# MF
clusterProfiler_MF_Stage_I   <- groupGO(gene     = df_expr_stage_I$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "MF", level    = 2,   readable = TRUE)@result
clusterProfiler_MF_Stage_II  <- groupGO(gene     = df_expr_stage_II$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "MF", level    = 2,   readable = TRUE)@result
clusterProfiler_MF_Stage_III <- groupGO(gene     = df_expr_stage_III$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "MF", level    = 2,   readable = TRUE)@result

# BP
clusterProfiler_BP_Stage_I   <- groupGO(gene     = df_expr_stage_I$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "BP", level    = 2,   readable = TRUE)@result
clusterProfiler_BP_Stage_II  <- groupGO(gene     = df_expr_stage_II$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "BP", level    = 2,   readable = TRUE)@result
clusterProfiler_BP_Stage_III <- groupGO(gene     = df_expr_stage_III$ENTREZID, OrgDb    = org.Hs.eg.db, ont      = "BP", level    = 2,   readable = TRUE)@result



ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
