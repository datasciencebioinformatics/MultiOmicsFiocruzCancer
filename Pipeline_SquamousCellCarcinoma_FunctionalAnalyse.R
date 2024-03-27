#######################################################################################################################################
# Path to files of selected_genes                                                                                                     # 
selected_genes_Stage_I_file       <-"/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_I.tsv"                            #
selected_genes_Stage_II_file      <-"/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_II.tsv"                           #
selected_genes_Stage_III_file     <-"/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_III.tsv"                          #
selected_genes_Stage_I_III_file   <-"/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_I_III.tsv"                        #
#######################################################################################################################################
# Load data                                                                                                                           #
selected_genes_Stage_I_data       <-read.table(file = selected_genes_Stage_I_file, sep = '\t', header = TRUE,fill=TRUE)               #
selected_genes_Stage_II_data      <-read.table(file = selected_genes_Stage_II_file, sep = '\t', header = TRUE,fill=TRUE)              #
selected_genes_Stage_III_data     <-read.table(file = selected_genes_Stage_III_file, sep = '\t', header = TRUE,fill=TRUE)             #
selected_genes_Stage_I_III_data   <-read.table(file = selected_genes_Stage_I_III_file, sep = '\t', header = TRUE,fill=TRUE)           #
#######################################################################################################################################
library(ggVennDiagram)                                                                                                                #################################################################################################################################################################
stages_I_II_III<-ggVennDiagram(list(Stage_I    =selected_genes_Stage_I_data$Gene,Stage_III  =selected_genes_Stage_III_data$Gene,Stage_II   =selected_genes_Stage_II_data$Gene), label_alpha = 0) + scale_fill_viridis() + theme_bw() + ggtitle("Stages I, II and III")                                #
stages_I_vs_III<-ggVennDiagram(list(Stages_I    =selected_genes_Stage_I_data$Gene, Stages_III    =selected_genes_Stage_III_data$Gene, Stages_II_vs_III    =selected_genes_Stage_I_III_data$Gene), label_alpha = 0) + scale_fill_viridis() + theme_bw() + ggtitle("Stages I, III and Stage I vs III")  #
#######################################################################################################################################################################################################################################################################################################





