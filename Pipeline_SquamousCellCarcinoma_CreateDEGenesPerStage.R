####################################################################################################################
# Filter down colData:
# First,  only stages I and III
# Second, only stages II and III
# Third,  only stages II and III
####################################################################################################################
# Aalysed pairs
df_stages<-data.frame(First=c("Stage I","Stage I","Stage II"),Second=c("Stage II","Stage III","Stage III"))
####################################################################################################################
# Run DESeq2
dds_stages <- DESeqDataSetFromMatrix(countData = unstranded_data, colData=colData[colnames(unstranded_data),], design = ~  age_at_index +  gender +stages  )

# Run DESeq2
dds_stages <- DESeq(dds_stages)

# Obtain differential expression numbers
resultsNames(dds_stages)

# Df s6tages I
df_stage_I_vs_III<-data.frame(results(dds_stages,name="stages_Stage.III_vs_Stage.I"))
####################################################################################################################
# Run varianceStabilizingTransformation
vst_stages_sub<-varianceStabilizingTransformation(dds_stages, blind = TRUE, fitType = "parametric")
####################################################################################################################
# First, stage O
# First, stageI
# Sort table by abs(log2FoldChange) and -log(padj)
down_regulated_df_stage_I_vs_III<-df_stage_I_vs_III[which(df_stage_I_vs_III$log2FoldChange>=0), ]
up_regulated_df_stage_I_vs_III<-df_stage_I_vs_III[which(df_stage_I_vs_III$log2FoldChange<0), ]
####################################################################################################################
# First, stageI
# Sort table by abs(log2FoldChange) and -log(padj)
Normal_Tumor_up_sort_df_stage_I_vs_III<- up_regulated_df_stage_I_vs_III[order(up_regulated_df_stage_I_vs_III$padj), ]
Normal_Tumor_down_sort_df_stage_I_vs_III<- down_regulated_df_stage_I_vs_III[order(down_regulated_df_stage_I_vs_III$padj), ]
Normal_Tumor_up_sort_df_stage_I_vs_III<- up_regulated_df_stage_I_vs_III[order(up_regulated_df_stage_I_vs_III$padj), ]
Normal_Tumor_down_sort_df_stage_I_vs_III<- down_regulated_df_stage_I_vs_III[order(down_regulated_df_stage_I_vs_III$padj), ]
####################################################################################################################
# Remove NA rows
# First, stageI
Normal_Tumor_up_sort_df_stage_I_vs_III<-na.omit(Normal_Tumor_up_sort_df_stage_I_vs_III)
Normal_Tumor_down_sort_df_stage_I_vs_III<-na.omit(Normal_Tumor_down_sort_df_stage_I_vs_III)
####################################################################################################################
# Field for top 10 percent of sorted sample
# First, stageI
Normal_Tumor_up_sort_df_stage_I_vs_III$Normal_Tumor_sort_10.0<-FALSE
Normal_Tumor_down_sort_df_stage_I_vs_III$Normal_Tumor_sort_10.0<-FALSE
####################################################################################################################
# First, stageI
# Field for top 10 percent of sorted sample
Normal_Tumor_up_sort_df_stage_I_vs_III[1:(dim(Normal_Tumor_up_sort_df_stage_I_vs_III)[1]*0.10),"Normal_Tumor_sort_10.0"]<-TRUE
Normal_Tumor_down_sort_df_stage_I_vs_III[1:(dim(Normal_Tumor_down_sort_df_stage_I_vs_III)[1]*0.10),"Normal_Tumor_sort_10.0"]<-TRUE
####################################################################################################################
Normal_Tumor_sort_stages<-rbind(Normal_Tumor_up_sort_df_stage_I_vs_III,Normal_Tumor_down_sort_df_stage_I_vs_III)
####################################################################################################################
# First, stageI
# "Unchanged"
Normal_Tumor_sort_stages$Expression<-0
####################################################################################################################
# Set expression up
# First, stageI
Normal_Tumor_sort_stages[intersect(which(Normal_Tumor_sort_stages$Normal_Tumor_sort_10.0), which(Normal_Tumor_sort_stages$log2FoldChange < 0)),"Expression"]<--1
Normal_Tumor_sort_stages[intersect(which(Normal_Tumor_sort_stages$Normal_Tumor_sort_10.0), which(Normal_Tumor_sort_stages$log2FoldChange >= 0)),"Expression"]<-1
####################################################################################################################
# First, stageI
Normal_Tumor_sort_stages$Categories<-""
Normal_Tumor_sort_stages[which(Normal_Tumor_sort_stages$Expression==0),"Categories"]<-"Uncategorized"
Normal_Tumor_sort_stages[which(Normal_Tumor_sort_stages$Expression==1),"Categories"]<-"Up-regulated"
Normal_Tumor_sort_stages[which(Normal_Tumor_sort_stages$Expression==-1),"Categories"]<-"Down-regulated"
####################################################################################################################
# Save TSV file with genes from Stage3
write_tsv(Normal_Tumor_sort_stages, "/home/felipe/Documentos/LungPortal/output/genes_StageIII.tsv")
####################################################################################################################
# Create volcano plot
p1 <- ggplot(Normal_Tumor_sort_stages, aes(log2FoldChange, -log(padj,2))) +  geom_point(size = 2/5) +  theme_bw()

# The thresholds
threshold_padj<-min(-log(Normal_Tumor_sort_stages[Normal_Tumor_sort_stages$Categories!="Uncategorized","padj"]))
threshold_log2fc_up<-min(Normal_Tumor_sort_stages[Normal_Tumor_sort_stages$Categories=="Up-regulated","log2FoldChange"])
threshold_log2fc_down<-max(Normal_Tumor_sort_stages[Normal_Tumor_sort_stages$Categories=="Down-regulated","log2FoldChange"])

# Adding color to differentially expressed genes (DEGs)
p2 <- ggplot(Normal_Tumor_sort_stages, aes(log2FoldChange, -log(padj),color = Categories)) + geom_point(size = 2/5,aes(color = Categories))  +
  xlab(expression("log2FoldChange")) + 
  ylab(expression("-log(padj)")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_bw() + ggtitle(paste("DE Genes Stage I vs. Stages III \n10.0% of top up-regulated, 10.0% of top up-regulated (padj)\n",paste("Up-regulated :",sum(Normal_Tumor_sort_stages$Categories=="Up-regulated"),"Down-regulated :",sum(Normal_Tumor_sort_stages$Categories=="Down-regulated"),sep=" "))) 

# Add treshold lines
p2 <- p2 + geom_hline(yintercept=threshold_padj ,linetype = 'dashed') + geom_vline(xintercept=threshold_log2fc_up ,linetype = 'dashed') + geom_vline(xintercept=threshold_log2fc_down ,linetype = 'dashed')

# Selected genes
selected_genes<-rownames(Normal_Tumor_sort_stages[Normal_Tumor_sort_stages$Categories!="Uncategorized",])
 
# Obtain differential expression numbers
pca_normal_stageIII<-plotPCA(vst_stages_sub[selected_genes,], intgroup="stages") + theme_bw() + ggtitle("DE Genes of Stage I vs. Stages III") + theme(legend.position='bottom')

Normal_Tumor_sort_stages<-Normal_Tumor_sort_stages[Normal_Tumor_sort_stages$Categories!="Uncategorized",]

# Change histogram plot fill colors by groups
padj_histogram_stages<-ggplot(Normal_Tumor_sort_stages, aes(x=-log(padj), fill=Categories, color=Categories)) +  geom_histogram(position="identity") + scale_fill_manual(values = c("dodgerblue3", "firebrick3"))  + theme_bw() 
#######################################################################################################################
# Remove samples from Stage III and plot again
pca_normal_stages_I_III<-plotPCA(vst_stages_sub[selected_genes,colData[colData$stages!="Stage II","patient_id"]], intgroup="stages") + theme_bw() + ggtitle("DE Genes of Stage I vs. Stage III") + theme(legend.position='bottom')
#######################################################################################################################
# FindClusters_resolution
png(filename=paste(output_dir,"Volcano_Plot_Normal_Tumor_Stage1_3.png",sep=""), width = 28, height = 24, res=600, units = "cm")
	pca_plots<-grid.arrange( p2, padj_histogram_stages,pca_normal_stageIII,pca_normal_stages_I_III,  ncol = 2)
dev.off()

# Save TSV file with genes from Stage1
# Save TSV file with genes from Stage1
write_tsv(Normal_Tumor_sort_stage_III[Normal_Tumor_sort_stage_III$Categories!="Uncategorized",], "/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_I_III.tsv")
########################################################################################################################
resultsNames(dds_stages)
# Df s6tages I
res<-results(dds_stages,name="stages_Stage.III_vs_Stage.I")
# FindClusters_resolution
png(filename=paste(output_dir,"Teste.png",sep=""), width = 28, height = 24, res=600, units = "cm")
	plotCounts(dds_stages, gene="ENSG00000169344.15", intgroup="stages")
dev.off()
