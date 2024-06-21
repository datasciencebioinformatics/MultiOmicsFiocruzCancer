##################################################################################################################################################
# Remove empty line
genes_stage_I<-genes_stage_I[genes_stage_I!=""]
genes_stage_II<-genes_stage_II[genes_stage_II!=""]
genes_stage_III<-genes_stage_III[genes_stage_III!=""]

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
# Take genes per each stage
genes_stage_I<-gsub(" ", "", unique(unique(strsplit(x=paste(unique(df_count_terms_selected$Genes_Stage_I),collapse=","),split=",",fixed=T))[[1]]))
genes_stage_II<-gsub(" ", "", unique(unique(strsplit(x=paste(unique(df_count_terms_selected$Genes_Stage_II),collapse=","),split=",",fixed=T))[[1]]))
genes_stage_III<-gsub(" ", "", unique(unique(strsplit(x=paste(unique(df_count_terms_selected$Genes_Stage_III),collapse=","),split=",",fixed=T))[[1]]))
genes_stage_I<-genes_stage_I[genes_stage_I!=""]
genes_stage_II<-genes_stage_II[genes_stage_II!=""]
genes_stage_III<-genes_stage_III[genes_stage_III!=""]
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
ids_translation      <-bitr(unstranded_data_samples$Gene, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL","ENSEMBL"), OrgDb="org.Hs.eg.db")
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
png(filename=paste(output_folder,"Plot_p_stage_I_paired.png",sep=""), width = 30, height = 30, res=600, units = "cm")
  grid.arrange(p_stage_I_paired,p_stage_II_paired,p_stage_III_paired)
dev.off()

# FindClusters_resolution
png(filename=paste(output_folder,"Plot_p_stage_I_unpaired.png",sep=""), width = 30, height = 30, res=600, units = "cm")
  grid.arrange(p_stage_I_unpaired,p_stage_II_unpaired,p_stage_III_unpaired)
dev.off()
