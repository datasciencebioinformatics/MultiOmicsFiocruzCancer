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
# Select only positive genes
selected_genes_Stage_I_pos<-selected_genes_Stage_I_data[which(selected_genes_Stage_I_data[,"log2FoldChange"]>0),]
selected_genes_Stage_II_pos<-selected_genes_Stage_II_data[which(selected_genes_Stage_II_data[,"log2FoldChange"]>0),]
selected_genes_Stage_III_pos<-selected_genes_Stage_III_data[which(selected_genes_Stage_III_data[,"log2FoldChange"]>0),]
selected_genes_Stage_I_III_pos<-selected_genes_Stage_I_III_data[which(selected_genes_Stage_I_III_data[,"log2FoldChange"]>0),]

# Calculate unique genes
unique_genes_Stage_I <-setdiff(selected_genes_Stage_I_pos$Gene,   c(selected_genes_Stage_II_pos$Gene,selected_genes_Stage_III_pos$Gene))
unique_genes_Stage_II <-setdiff(selected_genes_Stage_II_pos$Gene,  c(selected_genes_Stage_I_pos$Gene,selected_genes_Stage_III_pos$Gene))
unique_genes_Stage_III<-setdiff(selected_genes_Stage_III_pos$Gene, c(selected_genes_Stage_I_pos$Gene,selected_genes_Stage_II_pos$Gene))

# Set rownames
rownames(selected_genes_Stage_I_pos)<-selected_genes_Stage_I_pos$Gene
rownames(selected_genes_Stage_II_pos)<-selected_genes_Stage_II_pos$Gene
rownames(selected_genes_Stage_III_pos)<-selected_genes_Stage_III_pos$Gene
rownames(selected_genes_Stage_I_III_pos)<-selected_genes_Stage_I_III_pos$Gene
###################################################################################################################################################################################################################################################################################################
write_tsv(selected_genes_Stage_I_pos[unique_genes_Stage_I,],     "/home/felipe/Documentos/LungPortal/samples/pos_unique_genes_Stage_I.tsv")
write_tsv(selected_genes_Stage_II_pos[unique_genes_Stage_II,],   "/home/felipe/Documentos/LungPortal/samples/pos_unique_genes_Stage_II.tsv")
write_tsv(selected_genes_Stage_III_pos[unique_genes_Stage_III,],  "/home/felipe/Documentos/LungPortal/samples/pos_unique_genes_Stage_III.tsv")
write_tsv(selected_genes_Stage_I_III_pos, "/home/felipe/Documentos/LungPortal/samples/pos_unique_genes_Stages_I_III.tsv")
###################################################################################################################################################################################################################################################################################################
selected_genes_Stage_I_data       <- selected_genes_Stage_I_pos
selected_genes_Stage_II_data      <- selected_genes_Stage_II_pos
selected_genes_Stage_III_data     <- selected_genes_Stage_III_pos
selected_genes_Stage_I_III_data   <- selected_genes_Stage_I_III_pos

library(ggVennDiagram)
stages_I_II_III<-ggVennDiagram(list(Stage_I    =selected_genes_Stage_I_data$Gene,Stage_III  =selected_genes_Stage_III_data$Gene,Stage_II   =selected_genes_Stage_II_data$Gene), label_alpha = 0) + scale_fill_viridis() + theme_bw() + ggtitle("Stages I, II and III")                                #
stages_I_vs_III<-ggVennDiagram(list(Stages_I    =selected_genes_Stage_I_data$Gene, Stages_III    =selected_genes_Stage_III_data$Gene, Stages_I_vs_III    =selected_genes_Stage_I_III_data$Gene), label_alpha = 0) + scale_fill_viridis() + theme_bw() + ggtitle("Stages I, III and Stage I vs III")  #
#######################################################################################################################################################################################################################################################################################################
# Select stages 
stages_I_II_III<-ggVennDiagram(list(Stage_I    =selected_genes_Stage_I_data$Gene,Stage_III  =selected_genes_Stage_III_data$Gene,Stage_II   =selected_genes_Stage_II_data$Gene), label_alpha = 0) + scale_fill_viridis() + theme_bw() + ggtitle("Stages I, II and III")
stages_I_vs_III<-ggVennDiagram(list(Stages_I    =selected_genes_Stage_I_data$Gene, Stages_III    =selected_genes_Stage_III_data$Gene, Stages_II_vs_III    =selected_genes_Stage_I_III_data$Gene), label_alpha = 0) + scale_fill_viridis() + theme_bw() + ggtitle("Stages I, III and Stage I vs III")
###################################################################################################################################################################################################################################################################################################
# FindClusters_resolution
png(filename=paste(output_dir,"Veen_diagrams.png",sep=""), width = 14, height = 28, res=600, units = "cm")
  pca_plots<-grid.arrange( stages_I_II_III, stages_I_vs_III,  ncol = 1)
dev.off()




