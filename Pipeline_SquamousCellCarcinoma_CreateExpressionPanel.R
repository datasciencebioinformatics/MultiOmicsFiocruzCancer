#######################################################################################################################################
# Panel with gene expression per stage                                                                                                #
#######################################################################################################################################
# Path to files of selected_genes                                                                                                     # 
selected_genes_Stage_I_file       <-"/home/felipe/Documentos/LungPortal/output/uniq_selected_genes_Stage_pos_stage_I.tsv"             #
selected_genes_Stage_II_file      <-"/home/felipe/Documentos/LungPortal/output/uniq_selected_genes_Stage_pos_stage_II.tsv"            #
selected_genes_Stage_III_file     <-"/home/felipe/Documentos/LungPortal/output/uniq_selected_genes_Stage_pos_stage_III.tsv"           #
unstranded_file                   <-"/home/felipe/Documentos/LungPortal/samples/unstranded_data_id.tsv"                               #
df_log2foldchange_file            <-"/home/felipe/Documentos/LungPortal/samples/df_log2foldchange.tsv"                                #
#######################################################################################################################################
# Load data                                                                                                                           #
selected_genes_Stage_I_data       <-read.table(file = selected_genes_Stage_I_file, sep = '\t', header = TRUE,fill=TRUE)               #
selected_genes_Stage_II_data      <-read.table(file = selected_genes_Stage_II_file, sep = '\t', header = TRUE,fill=TRUE)              #
selected_genes_Stage_III_data     <-read.table(file = selected_genes_Stage_III_file, sep = '\t', header = TRUE,fill=TRUE)             #
unstranded_data                   <-read.table(file = unstranded_file, sep = '\t', header = TRUE,fill=TRUE)                           #
df_log2foldchange_data           <-read.table(file = df_log2foldchange_file, sep = '\t', header = TRUE,fill=TRUE)                     #
#######################################################################################################################################
# 3 genes specific to stage I in samples of Stage I, II and III                                                                       #
# 3 genes specific to stage II in samples of Stage I, II and III                                                                      #
# 3 genes specific to stage III in samples of Stage I, II and III                                                                     #
sel_genes_stage_I   <-selected_genes_Stage_I_data[order(-selected_genes_Stage_I_data$log2FoldChange),][1:3,]                          #
sel_genes_stage_II  <-selected_genes_Stage_II_data[order(-selected_genes_Stage_II_data$log2FoldChange),][1:3,]                        #
sel_genes_stage_III <-selected_genes_Stage_III_data[order(-selected_genes_Stage_III_data$log2FoldChange),][1:3,]                      #
#######################################################################################################################################
# Genes of each stage stored in colData                                                                                               #
sample_stage_I  <-colData[colData$stages=="Stage I","patient_id"]                                                                     #
sample_stage_II <-colData[colData$stages=="Stage II","patient_id"]                                                                    #
sample_stage_III<-colData[colData$stages=="Stage III","patient_id"]                                                                   #
#######################################################################################################################################
# Filter to keep only positive tumor/normal samples from unstranded_data.                                                             #
unstranded_data<-unstranded_data[avg_expression_pos$Gene,]                                                                            #
                                                                                                                                      #
# omit NA values from vector                                                                                                          #
unstranded_data <- na.omit(unstranded_data)                                                                                           #
                                                                                                                                      #
# Filter to keep only positive tumor/normal samples from colData_data.                                                                #
colData_data<-colData_data[colnames(unstranded_data),]                                                                                #
                                                                                                                                      #                                                                                                                                      #
# Create dds_stages with selected for the purpose of getting normalized counts                                                        #
dds_stages <- DESeqDataSetFromMatrix(countData = unstranded_data, colData=colData, design = ~  age_at_index +  gender +tissue_type  ) #
                                                                                                                                      #
# Estimate size factor                                                                                                                #
dds_stages <- estimateSizeFactors(dds_stages)                                                                                         #
                                                                                                                                      #
# Obtain normalized coutns                                                                                                            #
norm_counts<-counts(dds_stages, normalized = TRUE)                                                                                    #
#######################################################################################################################################
