library(ggVennDiagram)
library(viridis)

#######################################################################################################################################
# Path to files of selected_genes                                                                                                             # 
selected_genes_Stage_I_file       <-"/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_Stage_sample_stage_I.tsv"    #
selected_genes_Stage_II_file      <-"/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_Stage_sample_stage_II.tsv"   #
selected_genes_Stage_III_file     <-"/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_Stage_sample_stage_III.tsv"  #
#######################################################################################################################################
# Load data                                                                                                                           #
selected_genes_Stage_I_data       <-read.table(file = selected_genes_Stage_I_file, sep = '\t', header = TRUE,fill=TRUE)                        #
selected_genes_Stage_II_data      <-read.table(file = selected_genes_Stage_II_file, sep = '\t', header = TRUE,fill=TRUE)                       #
selected_genes_Stage_III_data     <-read.table(file = selected_genes_Stage_III_file, sep = '\t', header = TRUE,fill=TRUE)                      #
                                                                                                                                      #
# Set rownames                                                                                                                        #
rownames(selected_genes_Stage_I_data)<-selected_genes_Stage_I_data$gene                                                               #
rownames(selected_genes_Stage_II_data)<-selected_genes_Stage_II_data$gene                                                             #
rownames(selected_genes_Stage_III_data)<-selected_genes_Stage_III_data$gene                                                           #
#######################################################################################################################################
unique_stage_I=setdiff(selected_genes_Stage_I_data$gene, c(selected_genes_Stage_II_data$gene,selected_genes_Stage_III_data$gene))
unique_stage_II=setdiff(selected_genes_Stage_II_data$gene, c(selected_genes_Stage_I_data$gene,selected_genes_Stage_III_data$gene))
unique_stage_III=setdiff(selected_genes_Stage_III_data$gene, c(selected_genes_Stage_I_data$gene,selected_genes_Stage_II_data$gene))
#######################################################################################################################################
# Select stages 
stages_I_II_III_unique<-ggVennDiagram(list(Stage_I=unique_stage_I,Stage_II=unique_stage_II,Stage_III=unique_stage_III), label_alpha = 0) + scale_fill_viridis() + theme_bw() + ggtitle("Stages I, II and III")+ guides(fill="none")
stages_I_II_III<-ggVennDiagram(list(Stage_I=selected_genes_Stage_I_data$gene,Stage_II=selected_genes_Stage_II_data$gene,Stage_III=selected_genes_Stage_III_data$gene), label_alpha = 0) + scale_fill_viridis() + theme_bw() + ggtitle("Stages I, II and III, without overlap")+ guides(fill="none")

# FindClusters_resolution
png(filename=paste(output_dir,"Ven_Diagrams.png",sep=""), width = 28, height = 14, res=600, units = "cm")
  grid.arrange(stages_I_II_III, stages_I_II_III_unique, nrow = 1)
dev.off()
#############################################################################################################################
write_tsv(selected_genes_Stage_I_data[unique_stage_I,], "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_I.tsv")
write_tsv(selected_genes_Stage_II_data[unique_stage_II,], "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_II.tsv")
write_tsv(selected_genes_Stage_III_data[unique_stage_III,], "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedU_unique_III.tsv")
#############################################################################################################################