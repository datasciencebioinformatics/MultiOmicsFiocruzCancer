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
stages_I_II_III<-ggVennDiagram(list(Stage_I=selected_genes_Stage_I_data$gene,Stage_II=selected_genes_Stage_II_data$gene,Stage_III=selected_genes_Stage_III_data$gene), label_alpha = 0.9,set_color = c("grey50","grey50","grey50")) + scale_fill_gradient(low = "white", high = "white") + theme_bw() + ggtitle("Stage-specific genes")+ guides(fill="none")  + theme(panel.border = element_blank(), axis.text = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(), plot.background = element_rect(fill = "white"),  panel.background = element_rect(fill = "white"), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank())
stages_I_II_III<-ggVennDiagram(list("Stage I"=selected_genes_Stage_I_data$gene,"Stage II"=selected_genes_Stage_II_data$gene,"Stage III"=selected_genes_Stage_III_data$gene), label_alpha = 0.9,set_color = c("grey50","grey50","grey50")) + scale_fill_gradient(low = "white", high = "white") + theme_bw() + ggtitle("A")+ guides(fill="none")  + theme(panel.border = element_blank(), axis.text = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(), plot.background = element_rect(fill = "white"),  panel.background = element_rect(fill = "white"), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank())+labs( x ="", y = "")

# FindClusters_resolution
png(filename=paste(output_dir,"plot_ggVennDiagram_tumor_normal.png",sep=""), width = 16, height = 16, res=600, units = "cm")
  stages_I_II_III
dev.off()

#####################################################################################################################################################################################
# Tanspose RPKM table                                                                                                                                                               #
transporse_RPKM_table<-data.frame(t(unstranded_data_filter[,c(paired_sample_df$normal,paired_sample_df$tumor)]))                                                                                                                        #
                                                                                                                                                                                    #
# Calculate prcomp for stage                                                                                                                                                        #
pca_res_tumor_normal   <- prcomp(transporse_RPKM_table, scale. = TRUE)                                                                                                              #
                                                                                                                                                                                    #
# Plot PCA tumor versus normal                                                                                                                                                      #
plot_res_tumor_normal <- autoplot(pca_res_tumor_normal, data = colData[rownames(transporse_RPKM_table),], colour = 'tumor_normal')+ theme_bw()  + theme(legend.position="bottom") + ggtitle("B")                                                                      #

# FindClusters_resolution
png(filename=paste(output_dir,"plot_res_tumor_normal.png",sep=""), width = 16, height = 16, res=600, units = "cm")
  plot_res_tumor_normal
dev.off()

#####################################################################################################################################################################################
# FindClusters_resolution
png(filename=paste(output_dir,"plot_res_tumor_normal.png",sep=""), width = 20, height = 20, res=600, units = "cm")
  grid.arrange(plot_res_tumor_normal)
dev.off()
#####################################################################################################################################################################################
# ENSG00000120217 CD274/PD-L1
# ENSG00000171094 CD246/ALK 
# ENSG00000157764	BRAF
# ENSG00000133703 KRAS
# ENSG00000146648 EGFR
# ENSG00000047936	ROS1 
# ENSG00000141510	TP53
# ENSG00000121879	PIK3CA
# ENSG00000143924	EML4

# Tumor genes
# c("ENSG00000120217", "ENSG00000171094","ENSG00000157764", "ENSG00000133703", "ENSG00000146648", "ENSG00000047936", "ENSG00000141510", "ENSG00000121879", "ENSG00000143924")

# Stats of tumor genes
genes_unique_stages_filtered[which(genes_unique_stages_filtered$ENSEMBL %in% c("ENSG00000120217", "ENSG00000171094","ENSG00000157764", "ENSG00000133703", "ENSG00000146648", "ENSG00000047936", "ENSG00000141510", "ENSG00000121879", "ENSG00000143924")),]

# change box plot line colors by groups
p_stage_tumor_paired<-ggplot(unstranded_data_samples[unstranded_data_samples$ENSEMBL %in% c("ENSG00000133703", "ENSG00000141510", "ENSG00000146648", "ENSG00000143924"),], aes(x=tissue_type, y=RPKM, fill=tissue_type)) +   geom_boxplot()+ facet_wrap(~SYMBOL, ncol = 4, scales="free")+ theme_bw()

# FindClusters_resolution
png(filename=paste(output_dir,"plot_res_tumor_normal.png",sep=""), width = 16, height = 10, res=600, units = "cm")
  p_stage_tumor_paired
dev.off()
