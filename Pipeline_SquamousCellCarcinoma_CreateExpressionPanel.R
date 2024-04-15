#######################################################################################################################################

#######################################################################################################################################
# Path to files of selected_genes                                                                                                     # 
selected_genes_Stage_I_file       <-"/home/felipe/Documentos/LungPortal/output/uniq_selected_genes_Stage_pos_stage_I.tsv"             #
selected_genes_Stage_II_file      <-"/home/felipe/Documentos/LungPortal/output/uniq_selected_genes_Stage_pos_stage_II.tsv"            #
selected_genes_Stage_III_file     <-"/home/felipe/Documentos/LungPortal/output/uniq_selected_genes_Stage_pos_stage_III.tsv"           #
unstranded_file                   <-"/home/felipe/Documentos/LungPortal/samples/unstranded_data_id.tsv"                               #
#######################################################################################################################################
# Load data                                                                                                                           #
selected_genes_Stage_I_data       <-read.table(file = selected_genes_Stage_I_file, sep = '\t', header = TRUE,fill=TRUE)               #
selected_genes_Stage_II_data      <-read.table(file = selected_genes_Stage_II_file, sep = '\t', header = TRUE,fill=TRUE)              #
selected_genes_Stage_III_data     <-read.table(file = selected_genes_Stage_III_file, sep = '\t', header = TRUE,fill=TRUE)             #
unstranded_data                   <-read.table(file = unstranded_file, sep = '\t', header = TRUE,fill=TRUE)                           #
#######################################################################################################################################
