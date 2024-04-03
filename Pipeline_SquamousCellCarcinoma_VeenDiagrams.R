library(ggVennDiagram)
library(viridis)
#######################################################################################################################################
# Path to files of selected_genes                                                                                                     # 
selected_genes_Stage_I_file       <-"/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_pos_stage_I.tsv"                  #
selected_genes_Stage_II_file      <-"/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_pos_stage_II.tsv"                 #
selected_genes_Stage_III_file     <-"/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_pos_stage_III.tsv"                #
#######################################################################################################################################
# Load data                                                                                                                           #
selected_genes_Stage_I_data       <-read.table(file = selected_genes_Stage_I_file, sep = '\t', header = TRUE,fill=TRUE)               #
selected_genes_Stage_II_data      <-read.table(file = selected_genes_Stage_II_file, sep = '\t', header = TRUE,fill=TRUE)              #
selected_genes_Stage_III_data     <-read.table(file = selected_genes_Stage_III_file, sep = '\t', header = TRUE,fill=TRUE)             #

# Set rownames
rownames(selected_genes_Stage_I_pos)<-selected_genes_Stage_I_pos$Gene
rownames(selected_genes_Stage_II_pos)<-selected_genes_Stage_II_pos$Gene
rownames(selected_genes_Stage_III_pos)<-selected_genes_Stage_III_pos$Gene
#######################################################################################################################################
# Select stages 
stages_I_II_III<-ggVennDiagram(list(Stage_I    =selected_genes_Stage_I_data$Gene,Stage_II  =selected_genes_Stage_II_data$Gene,Stage_III  =selected_genes_Stage_III_data$Gene), label_alpha = 0) + scale_fill_viridis() + theme_bw() + ggtitle("Stages I, II and III")
#######################################################################################################################################
# Calculate unique genes
unique_genes_Stage_I <-setdiff(selected_genes_Stage_I_pos$Gene,   c(selected_genes_Stage_II_pos$Gene,selected_genes_Stage_III_pos$Gene))
unique_genes_Stage_II <-setdiff(selected_genes_Stage_II_pos$Gene,  c(selected_genes_Stage_I_pos$Gene,selected_genes_Stage_III_pos$Gene))
unique_genes_Stage_III<-setdiff(selected_genes_Stage_III_pos$Gene, c(selected_genes_Stage_I_pos$Gene,selected_genes_Stage_II_pos$Gene))
###################################################################################################################################################################################################################################################################################################
write_tsv(selected_genes_Stage_I_pos[unique_genes_Stage_I,],     "/home/felipe/Documentos/LungPortal/output/uniq_selected_genes_Stage_pos_stage_I.tsv")
write_tsv(selected_genes_Stage_II_pos[unique_genes_Stage_II,],   "/home/felipe/Documentos/LungPortal/output/uniq_selected_genes_Stage_pos_stage_II.tsv")
write_tsv(selected_genes_Stage_III_pos[unique_genes_Stage_III,],  "/home/felipe/Documentos/LungPortal/output/uniq_selected_genes_Stage_pos_stage_III.tsv")
###################################################################################################################################################################################################################################################################################################
selected_genes_Stage_I_data       <- unique_genes_Stage_I
selected_genes_Stage_II_data      <- unique_genes_Stage_II
selected_genes_Stage_III_data     <- unique_genes_Stage_III
#######################################################################################################################################################################################################################################################################################################
# Select stages 
stages_I_II_III_unique<-ggVennDiagram(list(Stage_I    =selected_genes_Stage_I_data,Stage_III  =selected_genes_Stage_III_data,Stage_II   =selected_genes_Stage_II_data), label_alpha = 0) + scale_fill_viridis() + theme_bw() + ggtitle("Unique genes Stages I, II and III")
###################################################################################################################################################################################################################################################################################################
# FindClusters_resolution
png(filename=paste(output_dir,"Veen_diagrams.png",sep=""), width = 28, height = 14, res=600, units = "cm")
  pca_plots<-grid.arrange( stages_I_II_III,stages_I_II_III_unique,  ncol = 2)
dev.off()
