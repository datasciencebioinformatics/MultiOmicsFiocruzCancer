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
sample_stage_I  <-colData[colData$sub_stages=="Stage I","patient_id"] 
sample_stage_IA  <-colData[colData$sub_stages=="Stage IA","patient_id"] 
sample_stage_IB  <-colData[colData$sub_stages=="Stage IB","patient_id"] 
sample_stage_II  <-colData[colData$sub_stages=="Stage II","patient_id"] 
sample_stage_IIA  <-colData[colData$sub_stages=="Stage IIA","patient_id"] 
sample_stage_IIB  <-colData[colData$sub_stages=="Stage IIB","patient_id"] 
sample_stage_III  <-colData[colData$sub_stages=="Stage III","patient_id"] 
sample_stage_IIIA  <-colData[colData$sub_stages=="Stage IIIA","patient_id"] 
sample_stage_IIIB  <-colData[colData$sub_stages=="Stage IIIB","patient_id"] 
#######################################################################################################################################
# Filter to keep only positive tumor/normal samples from unstranded_data.                                                             #
#unstranded_data<-unstranded_data[avg_expression_pos$Gene,]                                                                           #
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
# First, table for the gene selction is stored.
# genes from Stage I, pvalue <0.05
df_log2foldchange_data_stage_I_vs_II_III<-df_log2foldchange_data[df_log2foldchange_data$StageI_StagesII_III_pvalue<0.05,]
df_log2foldchange_data_stage_II_vs_I_III<-df_log2foldchange_data[df_log2foldchange_data$StageII_StagesI_III_pvalue<0.05,]
df_log2foldchange_data_stage_III_vs_I_II<-df_log2foldchange_data[df_log2foldchange_data$StageIII_StagesI_II_pvalue<0.05,]

# genes from Stage I, positive
df_log2foldchange_data_stage_I_vs_II_III<-df_log2foldchange_data_stage_I_vs_II_III[df_log2foldchange_data_stage_I_vs_II_III$StageI_StagesII_III_log2foldchange>0,]
df_log2foldchange_data_stage_II_vs_I_III<-df_log2foldchange_data_stage_II_vs_I_III[df_log2foldchange_data_stage_II_vs_I_III$StageII_StagesI_III_log2foldchange>0,]
df_log2foldchange_data_stage_III_vs_I_II<-df_log2foldchange_data_stage_III_vs_I_II[df_log2foldchange_data_stage_III_vs_I_II$StageIII_StagesI_II_log2foldchange>0,]

# Order by StageI_StagesII_III_log2foldchange
df_log2foldchange_data_stage_I_vs_II_III<-df_log2foldchange_data_stage_I_vs_II_III[order(-df_log2foldchange_data_stage_I_vs_II_III$StageI_StagesII_III_log2foldchange),c("Gene","StageI_StagesII_III_pvalue","StageI_StagesII_III_log2foldchange")][1:3,]
df_log2foldchange_data_stage_II_vs_I_III<-df_log2foldchange_data_stage_II_vs_I_III[order(-df_log2foldchange_data_stage_II_vs_I_III$StageII_StagesI_III_log2foldchange),c("Gene","StageII_StagesI_III_pvalue","StageII_StagesI_III_log2foldchange")][1:3,]
df_log2foldchange_data_stage_III_vs_I_II<-df_log2foldchange_data_stage_III_vs_I_II[order(-df_log2foldchange_data_stage_III_vs_I_II$StageIII_StagesI_II_log2foldchange),c("Gene","StageIII_StagesI_II_pvalue","StageI_StagesII_III_log2foldchange")][1:3,]

# A data frame for gene expression, samples and stage                                                                                                              #
# Set rownames                                                                                                                                                     #
# One gene specific to stage I in samples of Stage I, II and III                                                                                                   #
# Comb_genes_from_Intersections_vs_Samples_from_stage                                                                                                              #
sel_genes_stage_I_sample_stage_I        <-  data.frame(DE_genes_stage="Stage I",stages="Stage I", Expr=melt(norm_counts[df_log2foldchange_data_stage_I_vs_II_III$Gene,sample_stage_I]))       #
sel_genes_stage_I_sample_stage_IA       <-  data.frame(DE_genes_stage="Stage I",stages="Stage IA", Expr=melt(norm_counts[df_log2foldchange_data_stage_I_vs_II_III$Gene,sample_stage_IA]))     #
sel_genes_stage_I_sample_stage_IB       <-  data.frame(DE_genes_stage="Stage I",stages="Stage IB", Expr=melt(norm_counts[df_log2foldchange_data_stage_I_vs_II_III$Gene,sample_stage_IB]))     #
sel_genes_stage_I_sample_stage_II       <-  data.frame(DE_genes_stage="Stage I",stages="Stage II", Expr=melt(norm_counts[df_log2foldchange_data_stage_I_vs_II_III$Gene,sample_stage_II]))     #
sel_genes_stage_I_sample_stage_IIA      <-  data.frame(DE_genes_stage="Stage I",stages="Stage IIA", Expr=melt(norm_counts[df_log2foldchange_data_stage_I_vs_II_III$Gene,sample_stage_IIA]))   #
sel_genes_stage_I_sample_stage_IIB      <-  data.frame(DE_genes_stage="Stage I",stages="Stage IIB", Expr=melt(norm_counts[df_log2foldchange_data_stage_I_vs_II_III$Gene,sample_stage_IIB]))   #
sel_genes_stage_I_sample_stage_III      <-  data.frame(DE_genes_stage="Stage I",stages="Stage III", Expr=melt(norm_counts[df_log2foldchange_data_stage_I_vs_II_III$Gene,sample_stage_III]))   #
sel_genes_stage_I_sample_stage_IIIA     <-  data.frame(DE_genes_stage="Stage I",stages="Stage IIIA", Expr=melt(norm_counts[df_log2foldchange_data_stage_I_vs_II_III$Gene,sample_stage_IIIA])) #
sel_genes_stage_I_sample_stage_IIIB     <-  data.frame(DE_genes_stage="Stage I",stages="Stage IIIB", Expr=melt(norm_counts[df_log2foldchange_data_stage_I_vs_II_III$Gene,sample_stage_IIIB])) #

# A data frame for gene expression, samples and stage                                                                                                              #
# Set rownames                                                                                                                                                     #
# One gene specific to stage I in samples of Stage I, II and III                                                                                                   #
# Comb_genes_from_Intersections_vs_Samples_from_stage                                                                                                              #
sel_genes_stage_II_sample_stage_I        <-  data.frame(DE_genes_stage="Stage II",stages="Stage I", Expr=melt(norm_counts[df_log2foldchange_data_stage_II_vs_I_III$Gene,sample_stage_I]))       #
sel_genes_stage_II_sample_stage_IA       <-  data.frame(DE_genes_stage="Stage II",stages="Stage IA", Expr=melt(norm_counts[df_log2foldchange_data_stage_II_vs_I_III$Gene,sample_stage_IA]))     #
sel_genes_stage_II_sample_stage_IB       <-  data.frame(DE_genes_stage="Stage II",stages="Stage IB", Expr=melt(norm_counts[df_log2foldchange_data_stage_II_vs_I_III$Gene,sample_stage_IB]))     #
sel_genes_stage_II_sample_stage_II       <-  data.frame(DE_genes_stage="Stage II",stages="Stage II", Expr=melt(norm_counts[df_log2foldchange_data_stage_II_vs_I_III$Gene,sample_stage_II]))     #
sel_genes_stage_II_sample_stage_IIA      <-  data.frame(DE_genes_stage="Stage II",stages="Stage IIA", Expr=melt(norm_counts[df_log2foldchange_data_stage_II_vs_I_III$Gene,sample_stage_IIA]))   #
sel_genes_stage_II_sample_stage_IIB      <-  data.frame(DE_genes_stage="Stage II",stages="Stage IIB", Expr=melt(norm_counts[df_log2foldchange_data_stage_II_vs_I_III$Gene,sample_stage_IIB]))   #
sel_genes_stage_II_sample_stage_III      <-  data.frame(DE_genes_stage="Stage II",stages="Stage III", Expr=melt(norm_counts[df_log2foldchange_data_stage_II_vs_I_III$Gene,sample_stage_III]))   #
sel_genes_stage_II_sample_stage_IIIA     <-  data.frame(DE_genes_stage="Stage II",stages="Stage IIIA", Expr=melt(norm_counts[df_log2foldchange_data_stage_II_vs_I_III$Gene,sample_stage_IIIA])) #
sel_genes_stage_II_sample_stage_IIIB     <-  data.frame(DE_genes_stage="Stage II",stages="Stage IIIB", Expr=melt(norm_counts[df_log2foldchange_data_stage_II_vs_I_III$Gene,sample_stage_IIIB])) #

# A data frame for gene expression, samples and stage                                                                                                              #
# Set rownames                                                                                                                                                     #
# One gene specific to stage I in samples of Stage I, II and III                                                                                                   #
# Comb_genes_from_Intersections_vs_Samples_from_stage                                                                                                              #
sel_genes_stage_III_sample_stage_I        <-  data.frame(DE_genes_stage="Stage III",stages="Stage I", Expr=melt(norm_counts[df_log2foldchange_data_stage_III_vs_I_II$Gene,sample_stage_I]))       #
sel_genes_stage_III_sample_stage_IA       <-  data.frame(DE_genes_stage="Stage III",stages="Stage IA", Expr=melt(norm_counts[df_log2foldchange_data_stage_III_vs_I_II$Gene,sample_stage_IA]))     #
sel_genes_stage_III_sample_stage_IB       <-  data.frame(DE_genes_stage="Stage III",stages="Stage IB", Expr=melt(norm_counts[df_log2foldchange_data_stage_III_vs_I_II$Gene,sample_stage_IB]))     #
sel_genes_stage_III_sample_stage_II       <-  data.frame(DE_genes_stage="Stage III",stages="Stage II", Expr=melt(norm_counts[df_log2foldchange_data_stage_III_vs_I_II$Gene,sample_stage_II]))     #
sel_genes_stage_III_sample_stage_IIA      <-  data.frame(DE_genes_stage="Stage III",stages="Stage IIA", Expr=melt(norm_counts[df_log2foldchange_data_stage_III_vs_I_II$Gene,sample_stage_IIA]))   #
sel_genes_stage_III_sample_stage_IIB      <-  data.frame(DE_genes_stage="Stage III",stages="Stage IIB", Expr=melt(norm_counts[df_log2foldchange_data_stage_III_vs_I_II$Gene,sample_stage_IIB]))   #
sel_genes_stage_III_sample_stage_III      <-  data.frame(DE_genes_stage="Stage III",stages="Stage III", Expr=melt(norm_counts[df_log2foldchange_data_stage_III_vs_I_II$Gene,sample_stage_III]))   #
sel_genes_stage_III_sample_stage_IIIA     <-  data.frame(DE_genes_stage="Stage III",stages="Stage IIIA", Expr=melt(norm_counts[df_log2foldchange_data_stage_III_vs_I_II$Gene,sample_stage_IIIA])) #
sel_genes_stage_III_sample_stage_IIIB     <-  data.frame(DE_genes_stage="Stage III",stages="Stage IIIB", Expr=melt(norm_counts[df_log2foldchange_data_stage_III_vs_I_II$Gene,sample_stage_IIIB])) #

# Concatenate sel_genes_stageI
sel_genes_stageI  <-rbind(sel_genes_stage_I_sample_stage_I,sel_genes_stage_I_sample_stage_IA,sel_genes_stage_I_sample_stage_IB,sel_genes_stage_I_sample_stage_II,sel_genes_stage_I_sample_stage_IIA,sel_genes_stage_I_sample_stage_IIB,sel_genes_stage_I_sample_stage_III,sel_genes_stage_I_sample_stage_IIIA,sel_genes_stage_I_sample_stage_IIIB)
sel_genes_stageII <-rbind(sel_genes_stage_II_sample_stage_I,sel_genes_stage_II_sample_stage_IA,sel_genes_stage_II_sample_stage_IB,sel_genes_stage_II_sample_stage_II,sel_genes_stage_II_sample_stage_IIA,sel_genes_stage_II_sample_stage_IIB,sel_genes_stage_II_sample_stage_III,sel_genes_stage_II_sample_stage_IIIA,sel_genes_stage_II_sample_stage_IIIB)
sel_genes_stageIII<-rbind(sel_genes_stage_III_sample_stage_I,sel_genes_stage_III_sample_stage_IA,sel_genes_stage_III_sample_stage_IB,sel_genes_stage_III_sample_stage_II,sel_genes_stage_III_sample_stage_IIA,sel_genes_stage_III_sample_stage_IIB,sel_genes_stage_III_sample_stage_III,sel_genes_stage_III_sample_stage_IIIA,sel_genes_stage_III_sample_stage_IIIB)
                                                                                                                                                                   #
# Rename collumns
colnames(sel_genes_stageI)<-c("DE_genes_stage","stages","Gene","Patient","Expression")
colnames(sel_genes_stageII)<-c("DE_genes_stage","stages","Gene","Patient","Expression")
colnames(sel_genes_stageIII)<-c("DE_genes_stage","stages","Gene","Patient","Expression")

# Merge tables
sel_genes_stage<-rbind(sel_genes_stageI,sel_genes_stageII,sel_genes_stageIII)
####################################################################################################################################################################
#
# plot smean.sdl computes the mean plus or minus a constant times the standard deviation.                                                                                                                                                                                                                 #
p1 <- ggplot(sel_genes_stageI, aes(x=stages, y=Expression, fill=stages)) + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),geom="crossbar", width=0.5) +#
     facet_grid(~Gene, scales="free") +   
    theme(legend.position="none") + scale_fill_brewer(palette="Set1") + theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +                                                 #
    ggtitle("Expression of genes of Stage I") 


# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"Panel_genes_Substage1.png",sep=""), width = 16, height = 16, res=600, units = "cm")                                                                                                    #
  p1
dev.off() 


# plot                                                                                                                                                                                                                 #
p2 <- ggplot(sel_genes_stage, aes(x=stages, y=Expression, fill=stages)) +                                                                                                                                   #    
    theme(legend.position="none") + scale_fill_brewer(palette="Set1") + theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +                                                 #
    ggtitle("Expression of genes per stage") + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),geom="crossbar", width=0.5) + facet_wrap(~DE_genes_stage, scales = "free")


# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"Panel_genes_Substage1.png",sep=""), width = 16, height = 16, res=600, units = "cm")                                                                                                    #
  p2
dev.off() 
