#######################################################################################################################################
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CalculateShannonEntropyFromPairedUp.R")
#######################################################################################################################################
# Path to files of selected_genes                                                                                                     # 
selected_genes_Stage_I_file       <-paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_Stage_","sample_stage_I",".tsv",sep="")       #
selected_genes_Stage_II_file      <-paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_Stage_","sample_stage_II",".tsv",sep="")      #
selected_genes_Stage_III_file     <-paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_Stage_","sample_stage_III",".tsv",sep="")     #
#######################################################################################################################################
# Load data                                                                                                                           #
selected_genes_Stage_I_data       <-read.table(file = selected_genes_Stage_I_file, sep = '\t', header = TRUE,fill=TRUE)               #
selected_genes_Stage_II_data      <-read.table(file = selected_genes_Stage_II_file, sep = '\t', header = TRUE,fill=TRUE)              #
selected_genes_Stage_III_data     <-read.table(file = selected_genes_Stage_III_file, sep = '\t', header = TRUE,fill=TRUE)             #
                                                                                                                                      #
# Set rownames                                                                                                                        #
rownames(selected_genes_Stage_I_data)<-selected_genes_Stage_I_data$gene                                                               #
rownames(selected_genes_Stage_II_data)<-selected_genes_Stage_II_data$gene                                                             #
rownames(selected_genes_Stage_III_data)<-selected_genes_Stage_III_data$gene                                                           #
#####################################################################################################################################################################################
unique_stage_I  =intersect(setdiff(selected_genes_Stage_I_data$gene, c(selected_genes_Stage_II_data$gene,selected_genes_Stage_III_data$gene)),selected_genes_Stage_I_data$gene)     #
unique_stage_II =intersect(setdiff(selected_genes_Stage_II_data$gene, c(selected_genes_Stage_I_data$gene,selected_genes_Stage_III_data$gene)),selected_genes_Stage_II_data$gene)    #
unique_stage_III=intersect(setdiff(selected_genes_Stage_III_data$gene, c(selected_genes_Stage_I_data$gene,selected_genes_Stage_II_data$gene)),selected_genes_Stage_III_data$gene)   #
#####################################################################################################################################################################################
# Veen diagram
stages_I_II_III<-ggVennDiagram(list(Stage_I=selected_genes_Stage_I_data$gene,Stage_II=selected_genes_Stage_II_data$gene,Stage_III=selected_genes_Stage_III_data$gene), label_alpha = 0.9,set_color = c("grey50","grey50","grey50")) + scale_fill_gradient(low = "white", high = "white") + theme_bw() + ggtitle("Stages I, II and III")+ guides(fill="none")
#####################################################################################################################################################################################
# Tanspose RPKM table                                                                                                                                                               #
transporse_RPKM_table<-data.frame(t(unstranded_data_filter[,c(paired_sample_df$normal,paired_sample_df$tumor)]))                                                                                                                        #
                                                                                                                                                                                    #
# Calculate prcomp for stage                                                                                                                                                        #
pca_res_tumor_normal   <- prcomp(transporse_RPKM_table, scale. = TRUE)                                                                                                              #
                                                                                                                                                                                    #
# Plot PCA tumor versus normal                                                                                                                                                      #
plot_res_tumor_normal <- autoplot(pca_res_tumor_normal, data = colData[rownames(transporse_RPKM_table),], colour = 'tumor_normal')+ theme_bw()                                                                        #
#####################################################################################################################################################################################
# FindClusters_resolution
png(filename=paste(output_dir,"plot_res_tumor_normal.png",sep=""), width = 14, height = 14, res=600, units = "cm")
  plot_res_tumor_normal
dev.off()
#####################################################################################################################################################################################
# Selection of biomarkers
# Take genes of stage I
# Sort by vector name [z] then [x]
# Take genes of stage I
# Stage specific genes
stage_I_specific<-genes_unique_stages_stage_specific[genes_unique_stages_stage_specific[,"Stage_I.sig"]=="stage-specific",]
stage_II_specific<-genes_unique_stages_stage_specific[genes_unique_stages_stage_specific[,"Stage_II.sig"]=="stage-specific",]
stage_III_specific<-genes_unique_stages_stage_specific[genes_unique_stages_stage_specific[,"Stage_III.sig"]=="stage-specific",]

# Filter genes by FDR < 0.01
stage_I_specific<-stage_I_specific[stage_I_specific[,"Stage_I-normal.FDR"]<=0.01,]
stage_II_specific<-stage_II_specific[stage_II_specific[,"Stage_II-normal.FDR"]<=0.01,]
stage_III_specific<-stage_III_specific[stage_III_specific[,"Stage_III-normal.FDR"]<=0.01,]

# Filter genes by FDR < 0.01
stage_I_specific<-stage_I_specific[order(stage_I_specific[,"Stage_I-normal.log2fc"]),]
stage_II_specific<-stage_II_specific[order(stage_II_specific[,"Stage_II-normal.log2fc"]),]
stage_III_specific<-stage_III_specific[order(stage_III_specific[,"Stage_III-normal.log2fc"]),]
#####################################################################################################################################################################################
head(stage_I_specific,n=4)
head(stage_II_specific,n=4)
head(stage_III_specific,n=4)

# change box plot line colors by groups
p_stage_I_paired<-ggplot(unstranded_data_samples[unstranded_data_samples$ENSEMBL %in% head(stage_I_specific,n=4)[,c("ENSEMBL")],], aes(x=stages, y=RPKM, fill=tissue_type)) +   geom_boxplot()+ facet_wrap(~SYMBOL, ncol = 2, scales="free")+ theme_bw() + ggtitle("Stage I genes. Paired samples")
p_stage_I_unpaired<-ggplot(unstranded_data_samples_unapaired[unstranded_data_samples_unapaired$ENSEMBL %in% head(stage_I_specific,n=4)[,c("ENSEMBL")],], aes(x=stages, y=RPKM)) +   geom_boxplot()+ facet_wrap(~SYMBOL, ncol = 2, scales="free")+ theme_bw() + ggtitle("Stage I genes. Tumor samples") #  + stat_compare_means(comparisons = my_comparisons, method="t.test") 

# change box plot line colors by groups
p_stage_II_paired<-ggplot(unstranded_data_samples[unstranded_data_samples$ENSEMBL %in% head(stage_II_specific,n=4)[,c("ENSEMBL")],], aes(x=stages, y=RPKM, fill=tissue_type)) +   geom_boxplot()+ facet_wrap(~SYMBOL, ncol = 2, scales="free")+ theme_bw() + ggtitle("Stage II genes. Paired samples")
p_stage_II_unpaired<-ggplot(unstranded_data_samples_unapaired[unstranded_data_samples_unapaired$ENSEMBL %in% head(stage_II_specific,n=4)[,c("ENSEMBL")],], aes(x=stages, y=RPKM)) +   geom_boxplot()+ facet_wrap(~SYMBOL, ncol = 2, scales="free")+ theme_bw() + ggtitle("Stage II genes. Tumor samples") #  + stat_compare_means(comparisons = my_comparisons, method="t.test") 

# change box plot line colors by groups
p_stage_III_paired<-ggplot(unstranded_data_samples[unstranded_data_samples$ENSEMBL %in% head(stage_III_specific,n=4)[,c("ENSEMBL")],], aes(x=stages, y=RPKM, fill=tissue_type)) +   geom_boxplot()+ facet_wrap(~SYMBOL, ncol = 2, scales="free")+ theme_bw() + ggtitle("Stage III genes. Paired samples")
p_stage_III_unpaired<-ggplot(unstranded_data_samples_unapaired[unstranded_data_samples_unapaired$ENSEMBL %in% head(stage_III_specific,n=4)[,c("ENSEMBL")],], aes(x=stages, y=RPKM)) +   geom_boxplot()+ facet_wrap(~SYMBOL, ncol = 4, scales="free")+ theme_bw() + ggtitle("Stage III genes. Tumor samples") #  + stat_compare_means(comparisons = my_comparisons, method="t.test") 


