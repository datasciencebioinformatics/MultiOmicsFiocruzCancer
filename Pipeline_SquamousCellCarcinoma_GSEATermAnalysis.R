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
expr_stage_I    <-unstranded_rpkm[genes_Stage_I$gene,]
expr_stage_II   <-unstranded_rpkm[genes_Stage_II$gene,]
expr_stage_III  <-unstranded_rpkm[genes_Stage_III$gene,]

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
# Here I must check if if I use all the msigdbr databases or any in particulart
# Run GSEA 
# First, all msigdbr
pathways_C2_CP        <- msigdbr(species = "Homo sapiens", category = "C2", subcategory="CP")
pathways_C4_CGN       <- msigdbr(species = "Homo sapiens", category = "C4", subcategory="CGN")
pathways_C4_CM        <- msigdbr(species = "Homo sapiens", category = "C4", subcategory="CM")
pathways_C5_MF        <- msigdbr(species = "Homo sapiens", category = "C5", subcategory="MF")
pathways_C5_CC        <- msigdbr(species = "Homo sapiens", category = "C5", subcategory="CC")
pathways_C5_BP        <- msigdbr(species = "Homo sapiens", category = "C5", subcategory="BP")
pathways_C6           <- msigdbr(species = "Homo sapiens", category = "C6")
pathways_C7           <- msigdbr(species = "Homo sapiens", category = "C7")


# Split name of pathways
pathways_C2_CP <- split(as.character(pathways_C2_CP$entrez_gene), pathways_C2_CP$gs_name)
pathways_C4_CGN <- split(as.character(pathways_C4_CGN$entrez_gene), pathways_C4_CGN$gs_name)
pathways_C4_CM <- split(as.character(pathways_C4_CM$entrez_gene), pathways_C4_CM$gs_name)
pathways_C5_MF <- split(as.character(pathways_C5_MF$entrez_gene), pathways_C5_MF$gs_name)
pathways_C5_CC <- split(as.character(pathways_C5_CC$entrez_gene), pathways_C5_CC$gs_name)
pathways_C5_BP <- split(as.character(pathways_C5_BP$entrez_gene), pathways_C5_BP$gs_name)
pathways_C6 <- split(as.character(pathways_C7$entrez_gene), pathways_C6$gs_name)
pathways_C7 <- split(as.character(pathways_C7$entrez_gene), pathways_C7$gs_name)
#######################################################################################################################################
# Split name of pathways Stage I
stage_I_l2fc<-genes_Stage_I$log2change
names(stage_I_l2fc)<-df_expr_stage_I[genes_Stage_I$gene,"ENTREZID"]
stage_I_l2fc<-stage_I_l2fc[which(!is.na(names(stage_I_l2fc)))]

# Split name of pathways Stage II
stage_II_l2fc<-genes_Stage_II$log2change
names(stage_II_l2fc)<-df_expr_stage_II[genes_Stage_II$gene,"ENTREZID"]
stage_II_l2fc<-stage_II_l2fc[which(!is.na(names(stage_II_l2fc)))]

# Split name of pathways Stage III
stage_III_l2fc<-genes_Stage_III$log2change
names(stage_III_l2fc)<-df_expr_stage_III[genes_Stage_III$gene,"ENTREZID"]
stage_III_l2fc<-stage_III_l2fc[which(!is.na(names(stage_III_l2fc)))]
#######################################################################################################################################
# Calculate fgsea C2_CP
fgsaRes_pathways_C2_CP_stage_I     <- fgsea(pathways_C2_CP, stage_I_l2fc, minSize=5, maxSize=500)
fgsaRes_pathways_C2_CP_stage_II    <- fgsea(pathways_C2_CP, stage_II_l2fc, minSize=5, maxSize=500)
fgsaRes_pathways_C2_CP_stage_III   <- fgsea(pathways_C2_CP, stage_III_l2fc, minSize=5, maxSize=500)

# Calculate fgsea pathways_C4_CGN
fgsaRes_pathways_C4_CGN_stage_I     <- fgsea(pathways_C4_CGN, stage_I_l2fc, minSize=5, maxSize=500)
fgsaRes_pathways_C4_CGN_stage_II    <- fgsea(pathways_C4_CGN, stage_II_l2fc, minSize=5, maxSize=500)
fgsaRes_pathways_C4_CGN_stage_III   <- fgsea(pathways_C4_CGN, stage_III_l2fc, minSize=5, maxSize=500)

# Calculate fgsea pathways_C4_CGN
fgsaRes_pathways_C4_CM_stage_I     <- fgsea(pathways_C4_CM, stage_I_l2fc, minSize=5, maxSize=500)
fgsaRes_pathways_C4_CM_stage_II    <- fgsea(pathways_C4_CM, stage_II_l2fc, minSize=5, maxSize=500)
fgsaRes_pathways_C4_CM_stage_III   <- fgsea(pathways_C4_CM, stage_III_l2fc, minSize=5, maxSize=500)

# Calculate fgsea pathways_C4_CGN
fgsaRes_pathways_C5_MF_stage_I     <- fgsea(pathways_C5_MF, stage_I_l2fc, minSize=5, maxSize=500)
fgsaRes_pathways_C5_MF_stage_II    <- fgsea(pathways_C5_MF, stage_II_l2fc, minSize=5, maxSize=500)
fgsaRes_pathways_C5_MF_stage_III   <- fgsea(pathways_C5_MF, stage_III_l2fc, minSize=5, maxSize=500)

# Calculate fgsea pathways_C5_CC
fgsaRes_pathways_C5_CC_stage_I     <- fgsea(pathways_C5_CC, stage_I_l2fc, minSize=5, maxSize=500)
fgsaRes_pathways_C5_CC_stage_II    <- fgsea(pathways_C5_CC, stage_II_l2fc, minSize=5, maxSize=500)
fgsaRes_pathways_C5_CC_stage_III   <- fgsea(pathways_C5_CC, stage_III_l2fc, minSize=5, maxSize=500)

# Calculate fgsea pathways_C5_BP
fgsaRes_pathways_C5_BP_stage_I     <- fgsea(pathways_C5_BP, stage_I_l2fc, minSize=5, maxSize=500)
fgsaRes_pathways_C5_BP_stage_II    <- fgsea(pathways_C5_BP, stage_II_l2fc, minSize=5, maxSize=500)
fgsaRes_pathways_C5_BP_stage_III   <- fgsea(pathways_C5_BP, stage_III_l2fc, minSize=5, maxSize=500)

# Calculate fgsea pathways_C6
fgsaRes_pathways_C6_stage_I     <- fgsea(pathways_C6, stage_I_l2fc, minSize=5, maxSize=500)
fgsaRes_pathways_C6_stage_II    <- fgsea(pathways_C6, stage_II_l2fc, minSize=5, maxSize=500)
fgsaRes_pathways_C6_stage_III   <- fgsea(pathways_C6, stage_III_l2fc, minSize=5, maxSize=500)

# Calculate fgsea pathways_C7
fgsaRes_pathways_C7_stage_I     <- fgsea(pathways_C7, stage_I_l2fc, minSize=5, maxSize=500)
fgsaRes_pathways_C7_stage_II    <- fgsea(pathways_C7, stage_II_l2fc, minSize=5, maxSize=500)
fgsaRes_pathways_C7_stage_III   <- fgsea(pathways_C7, stage_III_l2fc, minSize=5, maxSize=500)
#######################################################################################################################################
# Calculate fgsea
# Calculate fgsea
fgsaRes_pathways_C2_CP_stage_I     <- fgsaRes_pathways_C2_CP_stage_I[which(fgsaRes_pathways_C2_CP_stage_I$padj<0.01),]
fgsaRes_pathways_C2_CP_stage_II    <- fgsaRes_pathways_C2_CP_stage_II[which(fgsaRes_pathways_C2_CP_stage_II$padj<0.01),]
fgsaRes_pathways_C2_CP_stage_III   <- fgsaRes_pathways_C2_CP_stage_III[which(fgsaRes_pathways_C2_CP_stage_III$padj<0.01),]

# Calculate fgsea pathways_C4_CGN
fgsaRes_pathways_C4_CGN_stage_I     <- fgsaRes_pathways_C4_CGN_stage_I[which(fgsaRes_pathways_C4_CGN_stage_I$padj<0.01),]
fgsaRes_pathways_C4_CGN_stage_II    <- fgsaRes_pathways_C4_CGN_stage_II[which(fgsaRes_pathways_C4_CGN_stage_II$padj<0.01),]
fgsaRes_pathways_C4_CGN_stage_III   <- fgsaRes_pathways_C4_CGN_stage_III[which(fgsaRes_pathways_C4_CGN_stage_III$padj<0.01),]

# Calculate fgsea pathways_C4_CM
fgsaRes_pathways_C4_CM_stage_I     <- fgsaRes_pathways_C4_CM_stage_I[which(fgsaRes_pathways_C4_CM_stage_I$padj<0.01),]
fgsaRes_pathways_C4_CM_stage_II    <- fgsaRes_pathways_C4_CM_stage_II[which(fgsaRes_pathways_C4_CM_stage_II$padj<0.01),]
fgsaRes_pathways_C4_CM_stage_III   <- fgsaRes_pathways_C4_CM_stage_III[which(fgsaRes_pathways_C4_CM_stage_III$padj<0.01),]

# Calculate fgsea pathways_C4_CGN
fgsaRes_pathways_C5_MF_stage_I     <- fgsaRes_pathways_C5_MF_stage_I[which(fgsaRes_pathways_C5_MF_stage_I$padj<0.01),]
fgsaRes_pathways_C5_MF_stage_II    <- fgsaRes_pathways_C5_MF_stage_II[which(fgsaRes_pathways_C5_MF_stage_II$padj<0.01),]
fgsaRes_pathways_C5_MF_stage_III   <- fgsaRes_pathways_C5_MF_stage_III[which(fgsaRes_pathways_C5_MF_stage_III$padj<0.01),]

# Calculate fgsea pathways_C5_CC
fgsaRes_pathways_C5_CC_stage_I     <-   fgsaRes_pathways_C5_CC_stage_I[which(fgsaRes_pathways_C5_CC_stage_I$padj<0.01),]
fgsaRes_pathways_C5_CC_stage_II    <-   fgsaRes_pathways_C5_CC_stage_II[which(fgsaRes_pathways_C5_CC_stage_II$padj<0.01),]
fgsaRes_pathways_C5_CC_stage_III   <-   fgsaRes_pathways_C5_CC_stage_III[which(fgsaRes_pathways_C5_CC_stage_III$padj<0.01),]

# Calculate fgsea pathways_C5_BP
fgsaRes_pathways_C5_BP_stage_I     <-   fgsaRes_pathways_C5_BP_stage_I[which(fgsaRes_pathways_C5_BP_stage_I$padj<0.01),]
fgsaRes_pathways_C5_BP_stage_II    <-   fgsaRes_pathways_C5_BP_stage_II[which(fgsaRes_pathways_C5_BP_stage_II$padj<0.01),]
fgsaRes_pathways_C5_BP_stage_III   <-   fgsaRes_pathways_C5_BP_stage_III[which(fgsaRes_pathways_C5_BP_stage_III$padj<0.01),]

# Calculate fgsea
fgsaRes_pathways_C6_stage_I     <-      fgsaRes_pathways_C6_stage_I[which(fgsaRes_pathways_C6_stage_I$padj<0.01),]
fgsaRes_pathways_C6_stage_II    <-      fgsaRes_pathways_C6_stage_II[which(fgsaRes_pathways_C6_stage_II$padj<0.01),]
fgsaRes_pathways_C6_stage_III   <-      fgsaRes_pathways_C6_stage_III[which(fgsaRes_pathways_C6_stage_III$padj<0.01),]

# Calculate fgsea
fgsaRes_pathways_C7_stage_I     <-      fgsaRes_pathways_C7_stage_I[which(fgsaRes_pathways_C7_stage_I$padj<0.01),]
fgsaRes_pathways_C7_stage_II    <-      fgsaRes_pathways_C7_stage_II[which(fgsaRes_pathways_C7_stage_II$padj<0.01),]
fgsaRes_pathways_C7_stage_III   <-      fgsaRes_pathways_C7_stage_III[which(fgsaRes_pathways_C7_stage_III$padj<0.01),]
#######################################################################################################################################
fgsaRes_pathways_C2_CP_stage_I$Stage="Stage I"
fgsaRes_pathways_C2_CP_stage_II$Stage="Stage II"
fgsaRes_pathways_C2_CP_stage_III$Stage="Stage III"

fgsaRes_pathways_C4_CGN_stage_I$Stage="Stage I"
fgsaRes_pathways_C4_CGN_stage_II$Stage="Stage II"
fgsaRes_pathways_C4_CGN_stage_III$Stage="Stage III"

fgsaRes_pathways_C4_CM_stage_I$Stage="Stage I"
fgsaRes_pathways_C4_CM_stage_II$Stage="Stage II"
fgsaRes_pathways_C4_CM_stage_III$Stage="Stage III"

fgsaRes_pathways_C5_MF_stage_I$Stage="Stage I"
fgsaRes_pathways_C5_MF_stage_II$Stage="Stage II"
fgsaRes_pathways_C5_MF_stage_III$Stage="Stage III"

fgsaRes_pathways_C5_BP_stage_I$Stage="Stage I"
fgsaRes_pathways_C5_BP_stage_II$Stage="Stage II"
fgsaRes_pathways_C5_BP_stage_III$Stage="Stage III"

fgsaRes_pathways_C5_CC_stage_I$Stage="Stage I"
fgsaRes_pathways_C5_CC_stage_II$Stage="Stage II"
fgsaRes_pathways_C5_CC_stage_III$Stage="Stage III"

fgsaRes_pathways_C6_stage_I$Stage="Stage I"
fgsaRes_pathways_C6_stage_II$Stage="Stage II"
fgsaRes_pathways_C6_stage_III$Stage="Stage III"

fgsaRes_pathways_C7_stage_I$Stage="Stage I"
fgsaRes_pathways_C7_stage_II$Stage="Stage II"
fgsaRes_pathways_C7_stage_III$Stage="Stage III"

fgsaRes_pathways_C2_CP_all_stage <-rbind(fgsaRes_pathways_C2_CP_stage_I,fgsaRes_pathways_C2_CP_stage_II,fgsaRes_pathways_C2_CP_stage_III)
fgsaRes_pathways_C4_CGN_all_stage<-rbind(fgsaRes_pathways_C4_CGN_stage_I,fgsaRes_pathways_C4_CGN_stage_II,fgsaRes_pathways_C4_CGN_stage_III)
fgsaRes_pathways_C4_CM_all_stage<-rbind(fgsaRes_pathways_C4_CM_stage_I,fgsaRes_pathways_C4_CM_stage_II,fgsaRes_pathways_C4_CM_stage_III)
fgsaRes_pathways_C5_MF_all_stage<-rbind(fgsaRes_pathways_C5_MF_stage_I,fgsaRes_pathways_C5_MF_stage_II,fgsaRes_pathways_C5_MF_stage_III)
fgsaRes_pathways_C5_BP_all_stage<-rbind(fgsaRes_pathways_C5_BP_stage_I,fgsaRes_pathways_C5_BP_stage_II,fgsaRes_pathways_C5_BP_stage_III)
fgsaRes_pathways_C5_CC_all_stage<-rbind(fgsaRes_pathways_C5_CC_stage_I,fgsaRes_pathways_C5_CC_stage_II,fgsaRes_pathways_C5_CC_stage_III)
fgsaRes_pathways_C6_all_stage<-rbind(fgsaRes_pathways_C6_stage_I,fgsaRes_pathways_C6_stage_II,fgsaRes_pathways_C6_stage_III)
fgsaRes_pathways_C7_all_stage<-rbind(fgsaRes_pathways_C7_stage_I,fgsaRes_pathways_C7_stage_II,fgsaRes_pathways_C7_stage_III)

# Col names
col_names<-c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "Stage")

write.xlsx(data.frame(fgsaRes_pathways_C2_CP_all_stage)[,col_names], "fgsaRes_pathways_C2_CP_all_stage", file=paste(output_dir,"/clusters/fgsaRes_pathways_all.xlsx",sep=""),append = FALSE) # where x is a data.frame with a Date column.
write.xlsx(data.frame(fgsaRes_pathways_C4_CGN_all_stage)[,col_names], "fgsaRes_pathways_C4_CGN_all_stage", file=paste(output_dir,"/clusters/fgsaRes_pathways_all.xlsx",sep=""),append = TRUE) # where x is a data.frame with a Date column.
write.xlsx(data.frame(fgsaRes_pathways_C4_CM_all_stage)[,col_names], "fgsaRes_pathways_C4_CM_all_stage", file=paste(output_dir,"/clusters/fgsaRes_pathways_all.xlsx",sep=""),append = TRUE) # where x is a data.frame with a Date column.
write.xlsx(data.frame(fgsaRes_pathways_C5_MF_all_stage)[,col_names], "fgsaRes_pathways_C5_MF_all_stage", file=paste(output_dir,"/clusters/fgsaRes_pathways_all.xlsx",sep=""),append = TRUE) # where x is a data.frame with a Date column.
write.xlsx(data.frame(fgsaRes_pathways_C5_BP_all_stage[,col_names], "fgsaRes_pathways_C5_BP_all_stage", file=paste(output_dir,"/clusters/fgsaRes_pathways_all.xlsx",sep=""),append = TRUE) # where x is a data.frame with a Date column.
write.xlsx(data.frame(fgsaRes_pathways_C5_CC_all_stage[,col_names], "fgsaRes_pathways_C5_CC_all_stage", file=paste(output_dir,"/clusters/fgsaRes_pathways_all.xlsx",sep=""),append = TRUE) # where x is a data.frame with a Date column.
write.xlsx(data.frame(fgsaRes_pathways_C6_all_stage)[,col_names], "fgsaRes_pathways_C6_all_stage", file=paste(output_dir,"/clusters/fgsaRes_pathways_all.xlsx",sep=""),append = TRUE) # where x is a data.frame with a Date column.
write.xlsx(data.frame(fgsaRes_pathways_C7_all_stage)[,col_names], "fgsaRes_pathways_C7_all_stage", file=paste(output_dir,"/clusters/fgsaRes_pathways_all.xlsx",sep=""),append = TRUE) # where x is a data.frame with a Date column.
####################################################################################################################################### 
plot_pathways_C2_CP<-ggplot(data.frame(fgsaRes_pathways_C2_CP_all_stage), aes(x=substr(pathway,1,25), y=-log(padj,10), label=pathway)) +geom_bar(stat='identity', aes(fill=padj), width=.5) + coord_flip() + facet_grid(cols = vars(Stage))+ theme_bw() + ggtitle("pathways_C2_CP")
plot_pathways_C4_CGN<-ggplot(data.frame(fgsaRes_pathways_C4_CGN_all_stage), aes(x=substr(pathway,1,25), y=-log(padj,10), label=pathway)) +geom_bar(stat='identity', aes(fill=padj), width=.5) + coord_flip() + facet_grid(cols = vars(Stage))+ theme_bw() + ggtitle("pathways_C4_CGN")
plot_pathways_C4_CM<-ggplot(data.frame(fgsaRes_pathways_C4_CM_all_stage), aes(x=substr(pathway,1,25), y=-log(padj,10), label=pathway)) +geom_bar(stat='identity', aes(fill=padj), width=.5) + coord_flip() + facet_grid(cols = vars(Stage))+ theme_bw() + ggtitle("pathways_C4_CM")
plot_pathways_C5_MF<-ggplot(data.frame(fgsaRes_pathways_C5_MF_all_stage), aes(x=substr(pathway,1,25), y=-log(padj,10), label=pathway)) +geom_bar(stat='identity', aes(fill=padj), width=.5) + coord_flip() + facet_grid(cols = vars(Stage))+ theme_bw() + ggtitle("pathways_C5_MF")
plot_pathways_C5_BP<-ggplot(data.frame(fgsaRes_pathways_C5_BP_all_stage), aes(x=substr(pathway,1,25), y=-log(padj,10), label=pathway)) +geom_bar(stat='identity', aes(fill=padj), width=.5) + coord_flip() + facet_grid(cols = vars(Stage))+ theme_bw() + ggtitle("pathways_C5_BP")
plot_pathways_C5_CC<-ggplot(data.frame(fgsaRes_pathways_C5_CC_all_stage), aes(x=substr(pathway,1,25), y=-log(padj,10), label=pathway)) +geom_bar(stat='identity', aes(fill=padj), width=.5) + coord_flip() + facet_grid(cols = vars(Stage))+ theme_bw() + ggtitle("pathways_C5_CC")
plot_pathways_C6<-ggplot(data.frame(fgsaRes_pathways_C6_all_stage), aes(x=substr(pathway,1,25), y=-log(padj,10), label=pathway)) +geom_bar(stat='identity', aes(fill=padj), width=.5) + coord_flip() + facet_grid(cols = vars(Stage))+ theme_bw() + ggtitle("pathways_C6")
plot_pathways_C7<-ggplot(data.frame(fgsaRes_pathways_C7_all_stage), aes(x=substr(pathway,1,25), y=-log(padj,10), label=pathway)) +geom_bar(stat='identity', aes(fill=padj), width=.5) + coord_flip() + facet_grid(cols = vars(Stage))+ theme_bw() + ggtitle("pathways_C7")
#######################################################################################################################################
# FindClusters_resolution
png(filename=paste(output_folder,"plot_pathways_C2_CP.png",sep=""), width = 23, height = 16, res=600, units = "cm")
	plot_pathways_C2_CP
dev.off()
# FindClusters_resolution
png(filename=paste(output_folder,"plot_pathways_C4_CGN.png",sep=""), width = 23, height = 16, res=600, units = "cm")
	plot_pathways_C4_CGN
dev.off()
# FindClusters_resolution
png(filename=paste(output_folder,"plot_pathways_C4_CM.png",sep=""), width = 23, height = 16, res=600, units = "cm")
	plot_pathways_C4_CM
dev.off()
# FindClusters_resolution
png(filename=paste(output_folder,"plot_pathways_C5_MF.png",sep=""), width = 23, height = 16, res=600, units = "cm")
	plot_pathways_C5_MF
dev.off()
# FindClusters_resolution
png(filename=paste(output_folder,"plot_pathways_C5_BP.png",sep=""), width = 23, height = 16, res=600, units = "cm")
	plot_pathways_C5_BP
dev.off()
# FindClusters_resolution
png(filename=paste(output_folder,"plot_pathways_C5_CC.png",sep=""), width = 23, height = 16, res=600, units = "cm")
	plot_pathways_C5_CC
dev.off()	   
# FindClusters_resolution
png(filename=paste(output_folder,"plot_pathways_C6.png",sep=""), width = 23, height = 16, res=600, units = "cm")
	plot_pathways_C6
dev.off()	   
# FindClusters_resolution
png(filename=paste(output_folder,"plot_pathways_C7.png",sep=""), width = 23, height = 16, res=600, units = "cm")
	plot_pathways_C7
dev.off()	   	   
#######################################################################################################################################
plotGseaTable(pathways_C2_CP[fgsaRes_pathways_C2_CP_stage_I$pathway], stage_I_l2fc,              fgsaRes_pathways_C2_CP_stage_I, gseaParam=0.5)  
#######################################################################################################################################
