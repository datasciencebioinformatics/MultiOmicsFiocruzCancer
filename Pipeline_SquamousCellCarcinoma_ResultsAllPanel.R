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
