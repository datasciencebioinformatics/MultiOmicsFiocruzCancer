##################################################################################################################################################
# Remove empty line
genes_stage_I<-annotation_stage_I$Symbol
genes_stage_II<-annotation_stage_II$Symbol
genes_stage_III<-annotation_stage_III$Symbol

# Select only tumor metadata
colDta_tumor<-colData[colData$tissue_type=="Tumor",]

# Select only tumor samples id's
tumor_samples<-colDta_tumor$patient_id

samples_normal<-colData[paired_sample_df$normal,]
samples_tumor<-colData[paired_sample_df$tumor,]

# Select samples of each stage stored in colData                                                                                             #
sample_stage_I_tumor  <-samples_tumor[samples_tumor$stages=="Stage I",c("tissue_type","patient_id","stages")]                                #
sample_stage_I_normal  <-samples_normal[samples_normal$stages=="Stage I",c("tissue_type","patient_id","stages")]                             #
sample_stage_II_tumor  <-samples_tumor[samples_tumor$stages=="Stage II",c("tissue_type","patient_id","stages")]                              #
sample_stage_II_normal  <-samples_normal[samples_normal$stages=="Stage II",c("tissue_type","patient_id","stages")]                           #
sample_stage_III_tumor  <-samples_tumor[samples_tumor$stages=="Stage III",c("tissue_type","patient_id","stages")]                            #
sample_stage_III_normal  <-samples_normal[samples_normal$stages=="Stage III",c("tissue_type","patient_id","stages")]                         ##############################
                                                                                                                                                                          #
# Select samples                                                                                                                                                          #
sample_stage_all_samples<-rbind(sample_stage_I_tumor,sample_stage_I_normal,sample_stage_II_tumor,sample_stage_II_normal,sample_stage_III_tumor,sample_stage_III_normal)   #
############################################################################################################################################################################
# ids_stage_I - all ENSEMBL anotated using bitr
ids_stage_I      <-bitr(genes_unique_Stage_I$gene_id, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
ids_stage_II     <-bitr(genes_unique_Stage_II$gene_id, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
ids_stage_III    <-bitr(genes_unique_Stage_III$gene_id, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
############################################################################################################################################################################
# Samples
unstranded_data_samples<-melt(t(unstranded_data_filter[,sample_stage_all_samples$patient_id]))
unstranded_data_samples_unapaired<-melt(t(unstranded_data_filter[,colDta_tumor$patient_id]))

# colnames(unstranded_data_samples)
colnames(unstranded_data_samples)<-c("patient_id","gene_id","RPKM")
colnames(unstranded_data_samples_unapaired)<-c("patient_id","gene_id","RPKM")

# unstranded_data_samples with sample_stage_all_samples
unstranded_data_samples<-merge(unstranded_data_samples,sample_stage_all_samples,by="patient_id")
unstranded_data_samples_unapaired<-merge(unstranded_data_samples_unapaired,sample_stage_all_samples,by="patient_id")

# colnames
colnames(unstranded_data_samples)<-c("patient_id","gene_id","RPKM","tissue_type","stages")
colnames(unstranded_data_samples_unapaired)<-c("patient_id","gene_id","RPKM","tissue_type","stages")

# genes_ids
df_id_conversion<-data.frame(Gene=c(),Gene_id=c())

# For each gene, add gene_id
for (gene_row in rownames(unstranded_data_filter))
{	      
    # Store gene id in the vector
    # Simply trim the gene id before the "." to save it in the ENSEML format
    rownames_id<-strsplit(gene_row, split = "\\.")[[1]][1]  
    df_id_conversion<-rbind(df_id_conversion,data.frame(Gene_id=gene_row,Gene=rownames_id))
}
# Set rownames
rownames(df_id_conversion)<-df_id_conversion$Gene_id

# Selected ids
selected_ids<-df_id_conversion[df_id_conversion$Gene %in% c(genes_unique_Stage_I$gene_id,genes_unique_Stage_II$gene_id,genes_unique_Stage_III$gene_id),]

# Take only selected entries
unstranded_data_samples<-unstranded_data_samples[which(unstranded_data_samples$gene_id %in% selected_ids$Gene_id),]
unstranded_data_samples_unapaired<-unstranded_data_samples_unapaired[which(unstranded_data_samples_unapaired$gene_id %in% selected_ids$Gene_id),]

# Set colnames
colnames(df_id_conversion)[1]<-c("gene_id")

# unstranded_data_samples
unstranded_data_samples<-merge(unstranded_data_samples,df_id_conversion,by="gene_id")
unstranded_data_samples_unapaired<-merge(unstranded_data_samples_unapaired,df_id_conversion,by="gene_id")
############################################################################################################################################################################
my_comparisons <- list( c("Stage I", "Stage II"), c("Stage I", "Stage III"), c("Stage II", "Stage III") )
############################################################################################################################################################################
# ids_stage_I - all ENSEMBL anotated using bitr
ids_translation               <-bitr(unstranded_data_samples$Gene, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL","ENSEMBL"), OrgDb="org.Hs.eg.db")
ids_translation_unpaired      <-bitr(unstranded_data_samples_unapaired$Gene, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL","ENSEMBL"), OrgDb="org.Hs.eg.db")

# unstranded_data_samples$Gene
unstranded_data_samples$ENSEMBL<-unstranded_data_samples$Gene
unstranded_data_samples_unapaired$ENSEMBL<-unstranded_data_samples_unapaired$Gene

# Merge symbols
unstranded_data_samples<-merge(unstranded_data_samples,ids_translation,by="ENSEMBL")
unstranded_data_samples_unapaired<-merge(unstranded_data_samples_unapaired,ids_translation,by="ENSEMBL")
###############################################################################################################
# Stage I
stage_I_selected_genes<-c("KPNA4","DHX36","CXADR")

# Stage II
stage_II_selected_genes<-c("QTRT1","MYBBP1A","THY1")

# Stage III
stage_III_selected_genes<-c("BRMS1","CPSF1","IDH1")
###############################################################################################################
# change box plot line colors by groups
p_stage_I_paired<-ggplot(unstranded_data_samples[unstranded_data_samples$SYMBOL %in% stage_I_selected_genes,], aes(x=stages, y=RPKM, fill=tissue_type)) +   geom_boxplot()+ facet_wrap(~SYMBOL, ncol = 3, scales="free")+ theme_bw() + ggtitle("Stage I genes. Paired samples")
p_stage_I_unpaired<-ggplot(unstranded_data_samples_unapaired[unstranded_data_samples_unapaired$SYMBOL %in% stage_I_selected_genes,], aes(x=stages, y=RPKM)) +   geom_boxplot()+ facet_wrap(~SYMBOL, ncol = 3, scales="free")+ theme_bw()  + stat_compare_means(comparisons = my_comparisons, method="t.test") + ggtitle("Stage I genes. Tumor samples")

# change box plot line colors by groups
p_stage_II_paired<-ggplot(unstranded_data_samples[unstranded_data_samples$SYMBOL %in% stage_II_selected_genes,], aes(x=stages, y=RPKM, fill=tissue_type)) +   geom_boxplot()+ facet_wrap(~SYMBOL, ncol = 3, scales="free")+ theme_bw() + ggtitle("Stage II genes. Paired samples")
p_stage_II_unpaired<-ggplot(unstranded_data_samples_unapaired[unstranded_data_samples_unapaired$SYMBOL %in% stage_II_selected_genes,], aes(x=stages, y=RPKM)) +   geom_boxplot()+ facet_wrap(~SYMBOL, ncol = 3, scales="free")+ theme_bw()  + stat_compare_means(comparisons = my_comparisons, method="t.test") + ggtitle("Stage II genes. Tumor samples")

# change box plot line colors by groups
p_stage_III_paired<-ggplot(unstranded_data_samples[unstranded_data_samples$SYMBOL %in% stage_III_selected_genes,], aes(x=stages, y=RPKM, fill=tissue_type)) +   geom_boxplot()+ facet_wrap(~SYMBOL, ncol = 3, scales="free")+ theme_bw() + ggtitle("Stage III genes. Paired samples")
p_stage_III_unpaired<-ggplot(unstranded_data_samples_unapaired[unstranded_data_samples_unapaired$SYMBOL %in% stage_III_selected_genes,], aes(x=stages, y=RPKM)) +   geom_boxplot()+ facet_wrap(~SYMBOL, ncol = 3, scales="free")+ theme_bw()  + stat_compare_means(comparisons = my_comparisons, method="t.test") + ggtitle("Stage III genes. Tumor samples")
############################################################################################################################################################################
# FindClusters_resolution
png(filename=paste(output_folder,"Plot_p_stage_selected_paired.png",sep=""), width = 30, height = 30, res=600, units = "cm")
  grid.arrange(p_stage_I_paired,p_stage_II_paired,p_stage_III_paired)
dev.off()

# FindClusters_resolution
png(filename=paste(output_folder,"Plot_p_stage_selected_unpaired.png",sep=""), width = 30, height = 30, res=600, units = "cm")
  grid.arrange(p_stage_I_unpaired,p_stage_II_unpaired,p_stage_III_unpaired)
dev.off()
############################################################################################################################################################################
# Log2foldchange (Tumor-Normal)
log2change_tumor_control_table<-log2change_tumor_control
log2change_tumor_control$Gene_id<-""

# Table with the paired
# log2change_tumor_control_paired<-log2change_tumor_control

# For each gene, add gene_id
for (gene_row in rownames(log2change_tumor_control_table))
{	      
  # Store gene id in the vector
  # Simply trim the gene id before the "." to save it in the ENSEML format
  rownames_id<-strsplit(gene_row, split = "\\.")[[1]][1]  
  
  # Addd gene id
  log2change_tumor_control_table[gene_row,"Gene_id"]<-rownames_id    
}
# Formart output tables
ids_all_tumor_sample      <-bitr(log2change_tumor_control_table$Gene_id, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")

# Duplicate field Gene_id
ids_all_tumor_sample$Gene_id<-ids_all_tumor_sample$ENSEMBL

# log2change_tumor_control
log2change_tumor_control_table<-merge(log2change_tumor_control_table,ids_all_tumor_sample,by="Gene_id", all = TRUE)
############################################################################################################################################################################
# Log2foldchange (Stage-Normal)
genes_unique_Stage_I     <-list_stage_specific_genes[["sample_stage_I"]]
genes_unique_Stage_II    <-list_stage_specific_genes[["sample_stage_II"]]
genes_unique_Stage_III   <-list_stage_specific_genes[["sample_stage_III"]]

genes_unique_Stage_I$ENSEMBL  <-""
genes_unique_Stage_II$ENSEMBL <-""
genes_unique_Stage_III$ENSEMBL<-""

genes_unique_Stage_I$Category<-"tumor-genes"
genes_unique_Stage_II$Category<-"tumor-genes"
genes_unique_Stage_III$Category<-"tumor-genes"

genes_unique_Stage_I[which(genes_unique_Stage_I$gene %in% unique_stage_I),"Category"]<-"stage-specific"
genes_unique_Stage_II[which(genes_unique_Stage_I$gene %in% unique_stage_II),"Category"]<-"stage-specific"
genes_unique_Stage_III[which(genes_unique_Stage_I$gene %in% unique_stage_III),"Category"]<-"stage-specific"

# Set colnames
colnames(log2change_tumor_control_paired)<-c("gene","tumor-normal.log2fc","Category","tumor-normal.pvalue","tumor-normal.FDR")

# Select collumns
log2change_tumor_control_selected<-log2change_tumor_control_paired[,c("gene","tumor-normal.log2fc","tumor-normal.pvalue","tumor-normal.FDR")]
############################################################################################################################################################################
# For each gene, add gene_id
for (gene_row in rownames(genes_unique_Stage_I))
{	      
  # Store gene id in the vector
  # Simply trim the gene id before the "." to save it in the ENSEML format
  rownames_id<-strsplit(gene_row, split = "\\.")[[1]][1]  
  
  # Addd gene id
  genes_unique_Stage_I[gene_row,"ENSEMBL"]<-rownames_id    
}
# For each gene, add gene_id
for (gene_row in rownames(genes_unique_Stage_II))
{	      
  # Store gene id in the vector
  # Simply trim the gene id before the "." to save it in the ENSEML format
  rownames_id<-strsplit(gene_row, split = "\\.")[[1]][1]  
  
  # Addd gene id
  genes_unique_Stage_II[gene_row,"ENSEMBL"]<-rownames_id    
}
# For each gene, add gene_id
for (gene_row in rownames(genes_unique_Stage_III))
{	      
  # Store gene id in the vector
  # Simply trim the gene id before the "." to save it in the ENSEML format
  rownames_id<-strsplit(gene_row, split = "\\.")[[1]][1]  
  
  # Addd gene id
  genes_unique_Stage_III[gene_row,"ENSEMBL"]<-rownames_id    
}
# For each gene, add gene_id
for (gene_row in rownames(log2change_tumor_control_selected))
{	      
  # Store gene id in the vector
  # Simply trim the gene id before the "." to save it in the ENSEML format
  rownames_id<-strsplit(gene_row, split = "\\.")[[1]][1]  
  
  # Addd gene id
  log2change_tumor_control_selected[gene_row,"ENSEMBL"]<-rownames_id    
}
############################################################################################################################################################################
# Log2foldchange (Stage-Normal)
colnames(genes_unique_Stage_I)   <-c("gene","Stage_I-normal.log2fc","Stage_I.sig","Stage_I-normal.pvalue","Stage_I-normal.FDR","ENSEMBL")
colnames(genes_unique_Stage_II)  <-c("gene","Stage_II-normal.log2fc","Stage_II.sig","Stage_II-normal.pvalue","Stage_II-normal.FDR","ENSEMBL")
colnames(genes_unique_Stage_III) <-c("gene","Stage_III-normal.log2fc","Stage_III.sig","Stage_III-normal.pvalue","Stage_III-normal.FDR","ENSEMBL")

# genes_unique_stages
genes_unique_stages  <-merge(merge(genes_unique_Stage_I,genes_unique_Stage_II,by="gene"),genes_unique_Stage_III,by="gene")

# Convert all symbols
ids_genes_unique_stages      <-bitr(genes_unique_stages$ENSEMBL, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL","ENSEMBL"), OrgDb="org.Hs.eg.db")

# genes_unique_stages
genes_unique_stages_filtered  <-merge(genes_unique_stages,ids_genes_unique_stages,by="ENSEMBL")

# genes_unique_stages
#log2change_tumor_control_paired <-merge(log2change_tumor_control_paired,ids_genes_unique_stages,by="ENSEMBL")

# genes_unique_stages
genes_unique_stages_filtered<-genes_unique_stages_filtered[,c("ENSEMBL","gene","ENTREZID","SYMBOL","Stage_I.sig","Stage_I-normal.log2fc","Stage_I-normal.pvalue","Stage_I-normal.FDR","Stage_II.sig","Stage_II-normal.log2fc","Stage_II-normal.pvalue","Stage_II-normal.FDR","Stage_III.sig","Stage_III-normal.log2fc","Stage_III-normal.pvalue","Stage_III-normal.FDR")]

# Merge tables
genes_unique_stages_complete<-merge(log2change_tumor_control_selected,genes_unique_stages_filtered,by="ENSEMBL")

# Select collumns
genes_unique_stages_filtered<-unique(genes_unique_stages_complete[,c("ENSEMBL", "ENTREZID", "SYMBOL" , "tumor-normal.log2fc","tumor-normal.pvalue" , "tumor-normal.FDR", "Stage_I.sig","Stage_I-normal.log2fc", "Stage_I-normal.pvalue", "Stage_I-normal.FDR", "Stage_II.sig", "Stage_II-normal.log2fc", "Stage_II-normal.pvalue", "Stage_II-normal.FDR", "Stage_III.sig", "Stage_III-normal.log2fc", "Stage_III-normal.pvalue", "Stage_III-normal.FDR")])

# Select stage-spefific
genes_unique_stages_stage_specific<-unique(genes_unique_stages_filtered[c(which(genes_unique_stages_filtered$Stage_I.sig=="stage-specific"),which(genes_unique_stages_filtered$Stage_II.sig=="stage-specific"),which(genes_unique_stages_filtered$Stage_III.sig=="stage-specific")),])
############################################################################################################################################################################
# log2fchange_results
write.xlsx(x=na.omit(genes_unique_stages_filtered), sheet="tumor genes",file=paste(output_dir,"log2fchange_results.xlsx",sep=""))
write.xlsx(x=na.omit(genes_unique_stages_stage_specific), sheet="stage-specific genes",file=paste(output_dir,"log2fchange_results.xlsx",sep=""),append = TRUE)
############################################################################################################################################################################


