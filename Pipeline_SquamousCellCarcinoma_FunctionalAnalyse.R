#######################################################################################################################################
# Path to files of selected_genes                                                                                                     # 
selected_genes_Stage_I_file       <-"/home/felipe/Documentos/LungPortal/samples/pos_unique_genes_Stage_I.tsv"                         #
selected_genes_Stage_II_file      <-"/home/felipe/Documentos/LungPortal/samples/pos_unique_genes_Stage_II.tsv"                        #
selected_genes_Stage_III_file     <-"/home/felipe/Documentos/LungPortal/samples/pos_unique_genes_Stage_III.tsv"                       #
#######################################################################################################################################
# Load data                                                                                                                           #
selected_genes_Stage_I_data       <-read.table(file = selected_genes_Stage_I_file, sep = '\t', header = TRUE,fill=TRUE)               #
selected_genes_Stage_II_data      <-read.table(file = selected_genes_Stage_II_file, sep = '\t', header = TRUE,fill=TRUE)              #
selected_genes_Stage_III_data     <-read.table(file = selected_genes_Stage_III_file, sep = '\t', header = TRUE,fill=TRUE)             #
#######################################################################################################################################
