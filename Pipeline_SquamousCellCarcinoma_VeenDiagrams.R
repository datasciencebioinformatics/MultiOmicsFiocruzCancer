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
rownames(selected_genes_Stage_I_data)<-selected_genes_Stage_I_data$Gene
rownames(selected_genes_Stage_II_data)<-selected_genes_Stage_II_data$Gene
rownames(selected_genes_Stage_III_data)<-selected_genes_Stage_III_data$Gene
#######################################################################################################################################
# Select stages 
stages_I_II_III<-ggVennDiagram(list(Stage_I    =selected_genes_Stage_I_data$Gene,Stage_II  =selected_genes_Stage_II_data$Gene,Stage_III  =selected_genes_Stage_III_data$Gene), label_alpha = 0) + scale_fill_viridis() + theme_bw() + ggtitle("Stages I, II and III")
#######################################################################################################################################
# Calculate unique genes
unique_genes_Stage_I <-setdiff(selected_genes_Stage_I_data$Gene,   c(selected_genes_Stage_II_data$Gene,selected_genes_Stage_III_data$Gene))
unique_genes_Stage_II <-setdiff(selected_genes_Stage_II_data$Gene,   c(selected_genes_Stage_I_data$Gene,selected_genes_Stage_III_data$Gene))
unique_genes_Stage_III<-setdiff(selected_genes_Stage_III_data$Gene, c(selected_genes_Stage_I_data$Gene,selected_genes_Stage_II_data$Gene))
###################################################################################################################################################################################################################################################################################################
write_tsv(selected_genes_Stage_I_data[unique_genes_Stage_I,],     "/home/felipe/Documentos/LungPortal/output/uniq_selected_genes_Stage_pos_stage_I.tsv")
write_tsv(selected_genes_Stage_II_data[unique_genes_Stage_II,],     "/home/felipe/Documentos/LungPortal/output/uniq_selected_genes_Stage_pos_stage_II.tsv")
write_tsv(selected_genes_Stage_III_data[unique_genes_Stage_III,],     "/home/felipe/Documentos/LungPortal/output/uniq_selected_genes_Stage_pos_stage_III.tsv")
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
###################################################################################################################################################################################################################################################################################################
# Load data                                                                                                                           #
selected_genes_Stage_I_data       <-read.table(file = selected_genes_Stage_I_file, sep = '\t', header = TRUE,fill=TRUE)               #
selected_genes_Stage_II_data      <-read.table(file = selected_genes_Stage_II_file, sep = '\t', header = TRUE,fill=TRUE)              #
selected_genes_Stage_III_data     <-read.table(file = selected_genes_Stage_III_file, sep = '\t', header = TRUE,fill=TRUE)             #

# Set rownames
rownames(selected_genes_Stage_I_data)<-selected_genes_Stage_I_data$Gene
rownames(selected_genes_Stage_II_data)<-selected_genes_Stage_II_data$Gene
rownames(selected_genes_Stage_III_data)<-selected_genes_Stage_III_data$Gene

# Calculate unique genes
intersect_I   <-intersect(intersect(selected_genes_Stage_I_data$Gene,selected_genes_Stage_I_data$Gene),selected_genes_Stage_III_data$Gene)
intersect_II  <-intersect(selected_genes_Stage_I_data$Gene,selected_genes_Stage_II_data$Gene)
intersect_III <-intersect(selected_genes_Stage_I_data$Gene,selected_genes_Stage_III_data$Gene)
intersect_IV  <-intersect(selected_genes_Stage_II_data$Gene,selected_genes_Stage_III_data$Gene)

write_tsv(data.frame(Gene=intersect_I),      "/home/felipe/Documentos/LungPortal/output/intersection_genes_pos_Stages_I_II_III.tsv")
write_tsv(data.frame(Gene=intersect_II),     "/home/felipe/Documentos/LungPortal/output/intersection_genes_pos_Stages_I_II.tsv")
write_tsv(data.frame(Gene=intersect_III),    "/home/felipe/Documentos/LungPortal/output/intersection_genes_pos_Stages_I_III.tsv")
write_tsv(data.frame(Gene=intersect_IV),     "/home/felipe/Documentos/LungPortal/output/intersection_genes_pos_Stages_II_III.tsv")
###################################################################################################################################################################################################################################################################################################
