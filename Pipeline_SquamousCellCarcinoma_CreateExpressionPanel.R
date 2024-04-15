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
#unstranded_data<-unstranded_data[avg_expression_pos$Gene,]                                                                            #
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
####################################################################################################################################################################
# genes from Stage I, pvalue <0.05
df_log2foldchange_data_stage_I_vs_II_III<-df_log2foldchange_data[df_log2foldchange_data$StageI_StagesII_III_pvalue<0.05,]

# genes from Stage I, positive
df_log2foldchange_data_stage_I_vs_II_III<-df_log2foldchange_data_stage_I_vs_II_III[df_log2foldchange_data_stage_I_vs_II_III$StageI_StagesII_III_log2foldchange>0,]

# Order by StageI_StagesII_III_log2foldchange
df_log2foldchange_data_stage_I_vs_II_III<-df_log2foldchange_data_stage_I_vs_II_III[order(-df_log2foldchange_data_stage_I_vs_II_III$StageI_StagesII_III_log2foldchange),c("Gene","StageI_StagesII_III_pvalue","StageI_StagesII_III_log2foldchange")][1:3,]


# A data frame for gene expression, samples and stage                                                                                                              #
# Set rownames                                                                                                                                                     #
# One gene specific to stage I in samples of Stage I, II and III                                                                                                   #
# Comb_genes_from_Intersections_vs_Samples_from_stage                                                                                                              #
sel_genes_stage_I_sample_stage_I      <-  data.frame(DE_genes_stage="Stage I",stages="Stage I", Expr=melt(norm_counts[df_log2foldchange_data_stage_I_vs_II_III$Gene,sample_stage_I]))     #
sel_genes_stage_I_sample_stage_II     <-  data.frame(DE_genes_stage="Stage I",stages="Stage II", Expr=melt(norm_counts[df_log2foldchange_data_stage_I_vs_II_III$Gene,sample_stage_II]))   #
sel_genes_stage_I_sample_stage_III    <-  data.frame(DE_genes_stage="Stage I",stages="Stage III", Expr=melt(norm_counts[df_log2foldchange_data_stage_I_vs_II_III$Gene,sample_stage_III])) #
                                                                                                                                                                   #
# Select genes from stage I                                                                                                                                        #
sel_genes_stageI<-rbind(sel_genes_stage_I_sample_stage_I,sel_genes_stage_I_sample_stage_II,sel_genes_stage_I_sample_stage_III)                                     #

# Rename collumns
colnames(sel_genes_stageI)<-c("DE_genes_stage","stages","Gene","Patient","Expression")
####################################################################################################################################################################
# plot                                                                                                                                                                                                                 #
p1 <- ggplot(sel_genes_stageI, aes(x=stages, y=Expression, fill=stages)) +                                                                                                                                   #
    geom_boxplot(varwidth = TRUE, outliers=FALSE) + facet_grid(~Gene) +                                                                                                                                                             #
    theme(legend.position="none") + geom_boxplot() + scale_fill_brewer(palette="Set1") + theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +                                                 #
    ggtitle("Expression of genes of Stage I") + stat_summary(fun=mean, colour="darkred", geom="point",  shape=18, size=3, show.legend=FALSE)

                                                                                                                                                                                                                            #
