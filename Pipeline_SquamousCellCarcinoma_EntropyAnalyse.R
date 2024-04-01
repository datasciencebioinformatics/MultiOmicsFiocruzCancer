#######################################################################################################################################
# To DO:
# Here 1 : try extremes of padj threshold range (0.001 - 0.05)
# Here 2 : try extremes of topNodes 10-100
# Here 3 : check if KS classiKS can be used as score instead of fisher qvalue
# Here 4 : change threshold parameter
# Important is to improve legibility of plot
library("DescTools")
#######################################################################################################################################
# Path to files of selected_genes                                                                                                     # 
selected_genes_Stage_I_file       <-"/home/felipe/Documentos/LungPortal/output/pos_unique_genes_Stage_I.tsv"                         #
selected_genes_Stage_II_file      <-"/home/felipe/Documentos/LungPortal/output/pos_unique_genes_Stage_II.tsv"                         #
selected_genes_Stage_III_file     <-"/home/felipe/Documentos/LungPortal/output/pos_unique_genes_Stage_III.tsv"                         #
#######################################################################################################################################
# Load data                                                                                                                           #
selected_genes_Stage_I_data       <-read.table(file = selected_genes_Stage_I_file, sep = '\t', header = TRUE,fill=TRUE)               #
selected_genes_Stage_II_data      <-read.table(file = selected_genes_Stage_II_file, sep = '\t', header = TRUE,fill=TRUE)              #
selected_genes_Stage_III_data     <-read.table(file = selected_genes_Stage_III_file, sep = '\t', header = TRUE,fill=TRUE)             #
#######################################################################################################################################
# Use function entropy
#######################################################################################################################################
