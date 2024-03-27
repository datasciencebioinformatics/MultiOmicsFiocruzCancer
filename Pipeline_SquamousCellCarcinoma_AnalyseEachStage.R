########################################################################################################################
# A panel to analyse differential expression comparing samples of each stage against all others stages.
########################################################################################################################
# Set colData
colData$stage_I   <- "Stages_II_III"
colData$stage_II  <- "Stages_I_III"
colData$stage_III <- "Stages_I_II"

# Each stage
colData$stage_I[which(colData$stages=="Stage I")]<-"Stage I"
colData$stage_II[which(colData$stages=="Stage II")]<-"Stage II"
colData$stage_III[which(colData$stages=="Stage III")]<-"Stage III"

# Removal of the normal samples and recalculate differential expression only with tumor samples
# Filter up only "Primary Tumor"
colData<-colData[colData$tumor_normal=="Primary Tumor",]

# Filter up only "Primary Tumor" patients
unstranded_data<-unstranded_data[,as.vector(colData$patient_id)]

# Run DESeq2
dds_stages_stage_I <- DESeqDataSetFromMatrix(countData = unstranded_data, colData=colData[colnames(unstranded_data),], design = ~  age_at_index +  gender +stage_I  )
dds_stages_stage_II <- DESeqDataSetFromMatrix(countData = unstranded_data, colData=colData[colnames(unstranded_data),], design = ~  age_at_index + gender +stage_II  )
dds_stages_stage_III <- DESeqDataSetFromMatrix(countData = unstranded_data, colData=colData[colnames(unstranded_data),], design = ~  age_at_index +  gender +stage_III  )

# Run DESeq2
dds_stages_stage_I <- DESeq(dds_stages_stage_I)
dds_stages_stage_II <- DESeq(dds_stages_stage_II)
dds_stages_stage_III <- DESeq(dds_stages_stage_III)

# Obtain differential expression numbers
resultsNames(dds_stages_stage_I)
resultsNames(dds_stages_stage_II)
resultsNames(dds_stages_stage_III)

# Df s6tages I
df_stage_I<-data.frame(results(dds_stages_stage_I,name="stage_I_Stages_II_III_vs_Stage.I"))
df_stage_II<-data.frame(results(dds_stages_stage_II,name="stage_II_Stages_I_III_vs_Stage.II"))
df_stage_III<-data.frame(results(dds_stages_stage_III,name="stage_III_Stages_I_II_vs_Stage.III"))

# Verify all components of PCA
# Run varianceStabilizingTransformation
vst_Stage_I_sub<-varianceStabilizingTransformation(dds_stages_stage_I, blind = TRUE, fitType = "parametric")
vst_Stage_II_sub<-varianceStabilizingTransformation(dds_stages_stage_II, blind = TRUE, fitType = "parametric")
vst_Stage_III_sub<-varianceStabilizingTransformation(dds_stages_stage_III, blind = TRUE, fitType = "parametric")
####################################################################################################################
# First, stage O
# First, stageI
# Sort table by abs(log2FoldChange) and -log(padj)
down_regulated_df_stage_I<-df_stage_I[which(df_stage_I$log2FoldChange>=0), ]
up_regulated_df_stage_I<-df_stage_I[which(df_stage_I$log2FoldChange<0), ]
####################################################################################################################
# First, stageI
# Sort table by abs(log2FoldChange) and -log(padj)
Normal_Tumor_up_sort_df_stage_I<- up_regulated_df_stage_I[order(up_regulated_df_stage_I$padj), ]
Normal_Tumor_down_sort_df_stage_I<- down_regulated_df_stage_I[order(down_regulated_df_stage_I$padj), ]
Normal_Tumor_up_sort_df_stage_I<- up_regulated_df_stage_I[order(up_regulated_df_stage_I$padj), ]
Normal_Tumor_down_sort_df_stage_I<- down_regulated_df_stage_I[order(down_regulated_df_stage_I$padj), ]
####################################################################################################################
# Remove NA rows
# First, stageI
Normal_Tumor_up_sort_df_stage_I<-na.omit(Normal_Tumor_up_sort_df_stage_I)
Normal_Tumor_down_sort_df_stage_I<-na.omit(Normal_Tumor_down_sort_df_stage_I)
####################################################################################################################
# Field for top 10 percent of sorted sample
# First, stageI
Normal_Tumor_up_sort_df_stage_I$Normal_Tumor_sort_10.0<-FALSE
Normal_Tumor_down_sort_df_stage_I$Normal_Tumor_sort_10.0<-FALSE
####################################################################################################################
# First, stageI
# Field for top 10 percent of sorted sample
Normal_Tumor_up_sort_df_stage_I[1:(dim(Normal_Tumor_up_sort_df_stage_I)[1]*0.10),"Normal_Tumor_sort_10.0"]<-TRUE
Normal_Tumor_down_sort_df_stage_I[1:(dim(Normal_Tumor_down_sort_df_stage_I)[1]*0.10),"Normal_Tumor_sort_10.0"]<-TRUE
####################################################################################################################
Normal_Tumor_sort_Stage_I<-rbind(Normal_Tumor_up_sort_df_stage_I,Normal_Tumor_down_sort_df_stage_I)
####################################################################################################################
# First, stageI
# "Unchanged"
Normal_Tumor_sort_Stage_I$Expression<-0
####################################################################################################################
# Set expression up
# First, stageI
Normal_Tumor_sort_Stage_I[intersect(which(Normal_Tumor_sort_Stage_I$Normal_Tumor_sort_10.0), which(Normal_Tumor_sort_Stage_I$log2FoldChange < 0)),"Expression"]<--1
Normal_Tumor_sort_Stage_I[intersect(which(Normal_Tumor_sort_Stage_I$Normal_Tumor_sort_10.0), which(Normal_Tumor_sort_Stage_I$log2FoldChange >= 0)),"Expression"]<-1

####################################################################################################################
# First, stageI
Normal_Tumor_sort_Stage_I$Categories<-""
Normal_Tumor_sort_Stage_I[which(Normal_Tumor_sort_Stage_I$Expression==0),"Categories"]<-"Uncategorized"
Normal_Tumor_sort_Stage_I[which(Normal_Tumor_sort_Stage_I$Expression==1),"Categories"]<-"Up-regulated"
Normal_Tumor_sort_Stage_I[which(Normal_Tumor_sort_Stage_I$Expression==-1),"Categories"]<-"Down-regulated"
####################################################################################################################
# Save TSV file with genes from Stage1
Normal_Tumor_sort_Stage_I$Gene<-rownames(Normal_Tumor_sort_Stage_I)
write_tsv(Normal_Tumor_sort_Stage_I, "/home/felipe/Documentos/LungPortal/output/genes_StageI.tsv")
####################################################################################################################
# Create volcano plot
p1 <- ggplot(Normal_Tumor_sort_Stage_I, aes(log2FoldChange, -log(padj,2))) +  geom_point(size = 2/5) +  theme_bw()

# The thresholds
threshold_padj<-min(-log(Normal_Tumor_sort_Stage_I[Normal_Tumor_sort_Stage_I$Categories!="Uncategorized","padj"]))
threshold_log2fc_up<-min(Normal_Tumor_sort_Stage_I[Normal_Tumor_sort_Stage_I$Categories=="Up-regulated","log2FoldChange"])
threshold_log2fc_down<-max(Normal_Tumor_sort_Stage_I[Normal_Tumor_sort_Stage_I$Categories=="Down-regulated","log2FoldChange"])

# Adding color to differentially expressed genes (DEGs)
p2 <- ggplot(Normal_Tumor_sort_Stage_I, aes(log2FoldChange, -log(padj),color = Categories)) + geom_point(size = 2/5,aes(color = Categories))  +
  xlab(expression("log2FoldChange")) + 
  ylab(expression("-log(padj)")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_bw() + ggtitle(paste("DE Genes Stage I vs. Stages II and III \n10.0% of top up-regulated, 10.0% of top up-regulated (padj)\n",paste("Up-regulated :",sum(Normal_Tumor_sort_Stage_I$Categories=="Up-regulated"),"Down-regulated :",sum(Normal_Tumor_sort_Stage_I$Categories=="Down-regulated"),sep=" "))) 

# Add treshold lines
p2 <- p2 + geom_hline(yintercept=threshold_padj ,linetype = 'dashed') + geom_vline(xintercept=threshold_log2fc_up ,linetype = 'dashed') + geom_vline(xintercept=threshold_log2fc_down ,linetype = 'dashed')

# Selected genes
selected_genes<-rownames(Normal_Tumor_sort_Stage_I[Normal_Tumor_sort_Stage_I$Categories!="Uncategorized",])
 
# Obtain differential expression numbers
pca_normal_stageI<-plotPCA(vst_Stage_I_sub[selected_genes,], intgroup="stage_I") + theme_bw() + ggtitle("DE Genes of Stage I") + theme(legend.position='bottom')

Normal_Tumor_sort_Stage_I<-Normal_Tumor_sort_Stage_I[Normal_Tumor_sort_Stage_I$Categories!="Uncategorized",]

# Change histogram plot fill colors by groups
padj_histogram_Stage_I<-ggplot(Normal_Tumor_sort_Stage_I, aes(x=-log(padj), fill=Categories, color=Categories)) +  geom_histogram(position="identity") + scale_fill_manual(values = c("dodgerblue3", "firebrick3"))  + theme_bw() 
#######################################################################################################################
# FindClusters_resolution
png(filename=paste(output_dir,"Volcano_Plot_Normal_Tumor_Stage1.png",sep=""), width = 28, height = 24, res=600, units = "cm")
	pca_plots<-grid.arrange( p2, padj_histogram_Stage_I,pca_normal_stageI,  ncol = 2)
dev.off()

# Save TSV file with genes from Stage1
write_tsv(Normal_Tumor_sort_Stage_I[Normal_Tumor_sort_Stage_I$Categories!="Uncategorized",], "/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_I.tsv")
########################################################################################################################
# First, stage O
# First, stageI
# Sort table by abs(log2FoldChange) and -log(padj)
down_regulated_df_stage_II<-df_stage_II[which(df_stage_II$log2FoldChange>=0), ]
up_regulated_df_stage_II<-df_stage_II[which(df_stage_II$log2FoldChange<0), ]
####################################################################################################################
# First, stageI
# Sort table by abs(log2FoldChange) and -log(padj)
Normal_Tumor_up_sort_df_stage_II<- up_regulated_df_stage_II[order(up_regulated_df_stage_II$padj), ]
Normal_Tumor_down_sort_df_stage_II<- down_regulated_df_stage_II[order(down_regulated_df_stage_II$padj), ]
Normal_Tumor_up_sort_df_stage_II<- up_regulated_df_stage_II[order(up_regulated_df_stage_II$padj), ]
Normal_Tumor_down_sort_df_stage_II<- down_regulated_df_stage_II[order(down_regulated_df_stage_II$padj), ]
####################################################################################################################
# Remove NA rows
# First, stageI
Normal_Tumor_up_sort_df_stage_II<-na.omit(Normal_Tumor_up_sort_df_stage_II)
Normal_Tumor_down_sort_df_stage_II<-na.omit(Normal_Tumor_down_sort_df_stage_II)
####################################################################################################################
# Field for top 10 percent of sorted sample
# First, stageI
Normal_Tumor_up_sort_df_stage_II$Normal_Tumor_sort_10.0<-FALSE
Normal_Tumor_down_sort_df_stage_II$Normal_Tumor_sort_10.0<-FALSE
####################################################################################################################
# First, stageI
# Field for top 10 percent of sorted sample
Normal_Tumor_up_sort_df_stage_II[1:(dim(Normal_Tumor_up_sort_df_stage_II)[1]*0.10),"Normal_Tumor_sort_10.0"]<-TRUE
Normal_Tumor_down_sort_df_stage_II[1:(dim(Normal_Tumor_down_sort_df_stage_II)[1]*0.10),"Normal_Tumor_sort_10.0"]<-TRUE
####################################################################################################################
Normal_Tumor_sort_stage_II<-rbind(Normal_Tumor_up_sort_df_stage_II,Normal_Tumor_down_sort_df_stage_II)
####################################################################################################################
# First, stageI
# "Unchanged"
Normal_Tumor_sort_stage_II$Expression<-0
####################################################################################################################
# Set expression up
# First, stageI
Normal_Tumor_sort_stage_II[intersect(which(Normal_Tumor_sort_stage_II$Normal_Tumor_sort_10.0), which(Normal_Tumor_sort_stage_II$log2FoldChange < 0)),"Expression"]<--1
Normal_Tumor_sort_stage_II[intersect(which(Normal_Tumor_sort_stage_II$Normal_Tumor_sort_10.0), which(Normal_Tumor_sort_stage_II$log2FoldChange >= 0)),"Expression"]<-1

####################################################################################################################
# First, stageI
Normal_Tumor_sort_stage_II$Categories<-""
Normal_Tumor_sort_stage_II[which(Normal_Tumor_sort_stage_II$Expression==0),"Categories"]<-"Uncategorized"
Normal_Tumor_sort_stage_II[which(Normal_Tumor_sort_stage_II$Expression==1),"Categories"]<-"Up-regulated"
Normal_Tumor_sort_stage_II[which(Normal_Tumor_sort_stage_II$Expression==-1),"Categories"]<-"Down-regulated"
####################################################################################################################
# Save TSV file with genes from Stage1
Normal_Tumor_sort_stage_II$Gene<-rownames(Normal_Tumor_sort_stage_II)
write_tsv(Normal_Tumor_sort_stage_II, "/home/felipe/Documentos/LungPortal/output/genes_StageII.tsv")
####################################################################################################################
# Create volcano plot
p1 <- ggplot(Normal_Tumor_sort_stage_II, aes(log2FoldChange, -log(padj,2))) +  geom_point(size = 2/5) +  theme_bw()

# The thresholds
threshold_padj<-min(-log(Normal_Tumor_sort_stage_II[Normal_Tumor_sort_stage_II$Categories!="Uncategorized","padj"]))
threshold_log2fc_up<-min(Normal_Tumor_sort_stage_II[Normal_Tumor_sort_stage_II$Categories=="Up-regulated","log2FoldChange"])
threshold_log2fc_down<-max(Normal_Tumor_sort_stage_II[Normal_Tumor_sort_stage_II$Categories=="Down-regulated","log2FoldChange"])

# Adding color to differentially expressed genes (DEGs)
p2 <- ggplot(Normal_Tumor_sort_stage_II, aes(log2FoldChange, -log(padj),color = Categories)) + geom_point(size = 2/5,aes(color = Categories))  +
  xlab(expression("log2FoldChange")) + 
  ylab(expression("-log(padj)")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_bw() + ggtitle(paste("DE Genes Stage II vs. Stages I and III \n10.0% of top up-regulated, 10.0% of top up-regulated (padj)\n",paste("Up-regulated :",sum(Normal_Tumor_sort_stage_II$Categories=="Up-regulated"),"Down-regulated :",sum(Normal_Tumor_sort_stage_II$Categories=="Down-regulated"),sep=" "))) 

# Add treshold lines
p2 <- p2 + geom_hline(yintercept=threshold_padj ,linetype = 'dashed') + geom_vline(xintercept=threshold_log2fc_up ,linetype = 'dashed') + geom_vline(xintercept=threshold_log2fc_down ,linetype = 'dashed')

# Selected genes
selected_genes<-rownames(Normal_Tumor_sort_stage_II[Normal_Tumor_sort_stage_II$Categories!="Uncategorized",])
 
# Obtain differential expression numbers
pca_normal_stageII<-plotPCA(vst_Stage_II_sub[selected_genes,], intgroup="stage_II") + theme_bw() + ggtitle("DE Genes of Stage II") + theme(legend.position='bottom')

Normal_Tumor_sort_stage_II<-Normal_Tumor_sort_stage_II[Normal_Tumor_sort_stage_II$Categories!="Uncategorized",]

# Change histogram plot fill colors by groups
padj_histogram_stage_II<-ggplot(Normal_Tumor_sort_stage_II, aes(x=-log(padj), fill=Categories, color=Categories)) +  geom_histogram(position="identity") + scale_fill_manual(values = c("dodgerblue3", "firebrick3"))  + theme_bw() 
#######################################################################################################################
# FindClusters_resolution
png(filename=paste(output_dir,"Volcano_Plot_Normal_Tumor_Stage2.png",sep=""), width = 28, height = 24, res=600, units = "cm")
	pca_plots<-grid.arrange( p2, padj_histogram_stage_II,pca_normal_stageII,  ncol = 2)
dev.off()

# Save TSV file with genes from Stage1
write_tsv(Normal_Tumor_sort_stage_II[Normal_Tumor_sort_stage_II$Categories!="Uncategorized",], "/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_II.tsv")
####################################################################################################################
# First, stage O
# First, stageI
# Sort table by abs(log2FoldChange) and -log(padj)
down_regulated_df_stage_III<-df_stage_III[which(df_stage_III$log2FoldChange>=0), ]
up_regulated_df_stage_III<-df_stage_III[which(df_stage_III$log2FoldChange<0), ]
####################################################################################################################
# First, stageI
# Sort table by abs(log2FoldChange) and -log(padj)
Normal_Tumor_up_sort_df_stage_III<- up_regulated_df_stage_III[order(up_regulated_df_stage_III$padj), ]
Normal_Tumor_down_sort_df_stage_III<- down_regulated_df_stage_III[order(down_regulated_df_stage_III$padj), ]
Normal_Tumor_up_sort_df_stage_III<- up_regulated_df_stage_III[order(up_regulated_df_stage_III$padj), ]
Normal_Tumor_down_sort_df_stage_III<- down_regulated_df_stage_III[order(down_regulated_df_stage_III$padj), ]
####################################################################################################################
# Remove NA rows
# First, stageI
Normal_Tumor_up_sort_df_stage_III<-na.omit(Normal_Tumor_up_sort_df_stage_III)
Normal_Tumor_down_sort_df_stage_III<-na.omit(Normal_Tumor_down_sort_df_stage_III)
####################################################################################################################
# Field for top 10 percent of sorted sample
# First, stageI
Normal_Tumor_up_sort_df_stage_III$Normal_Tumor_sort_10.0<-FALSE
Normal_Tumor_down_sort_df_stage_III$Normal_Tumor_sort_10.0<-FALSE
####################################################################################################################
# First, stageI
# Field for top 10 percent of sorted sample
Normal_Tumor_up_sort_df_stage_III[1:(dim(Normal_Tumor_up_sort_df_stage_III)[1]*0.10),"Normal_Tumor_sort_10.0"]<-TRUE
Normal_Tumor_down_sort_df_stage_III[1:(dim(Normal_Tumor_down_sort_df_stage_III)[1]*0.10),"Normal_Tumor_sort_10.0"]<-TRUE
####################################################################################################################
Normal_Tumor_sort_stage_III<-rbind(Normal_Tumor_up_sort_df_stage_III,Normal_Tumor_down_sort_df_stage_III)
####################################################################################################################
# First, stageI
# "Unchanged"
Normal_Tumor_sort_stage_III$Expression<-0
####################################################################################################################
# Set expression up
# First, stageI
Normal_Tumor_sort_stage_III[intersect(which(Normal_Tumor_sort_stage_III$Normal_Tumor_sort_10.0), which(Normal_Tumor_sort_stage_III$log2FoldChange < 0)),"Expression"]<--1
Normal_Tumor_sort_stage_III[intersect(which(Normal_Tumor_sort_stage_III$Normal_Tumor_sort_10.0), which(Normal_Tumor_sort_stage_III$log2FoldChange >= 0)),"Expression"]<-1

####################################################################################################################
# First, stageI
Normal_Tumor_sort_stage_III$Categories<-""
Normal_Tumor_sort_stage_III[which(Normal_Tumor_sort_stage_III$Expression==0),"Categories"]<-"Uncategorized"
Normal_Tumor_sort_stage_III[which(Normal_Tumor_sort_stage_III$Expression==1),"Categories"]<-"Up-regulated"
Normal_Tumor_sort_stage_III[which(Normal_Tumor_sort_stage_III$Expression==-1),"Categories"]<-"Down-regulated"
####################################################################################################################
# Save TSV file with genes from Stage3
Normal_Tumor_sort_stage_III$Gene<-rownames(Normal_Tumor_sort_stage_III)
write_tsv(Normal_Tumor_sort_stage_III, "/home/felipe/Documentos/LungPortal/output/genes_StageIII.tsv")
####################################################################################################################
# Create volcano plot
p1 <- ggplot(Normal_Tumor_sort_stage_III, aes(log2FoldChange, -log(padj,2))) +  geom_point(size = 2/5) +  theme_bw()

# The thresholds
threshold_padj<-min(-log(Normal_Tumor_sort_stage_III[Normal_Tumor_sort_stage_III$Categories!="Uncategorized","padj"]))
threshold_log2fc_up<-min(Normal_Tumor_sort_stage_III[Normal_Tumor_sort_stage_III$Categories=="Up-regulated","log2FoldChange"])
threshold_log2fc_down<-max(Normal_Tumor_sort_stage_III[Normal_Tumor_sort_stage_III$Categories=="Down-regulated","log2FoldChange"])

# Adding color to differentially expressed genes (DEGs)
p2 <- ggplot(Normal_Tumor_sort_stage_III, aes(log2FoldChange, -log(padj),color = Categories)) + geom_point(size = 2/5,aes(color = Categories))  +
  xlab(expression("log2FoldChange")) + 
  ylab(expression("-log(padj)")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_bw() + ggtitle(paste("DE Genes Stage III vs. Stages I and II \n10.0% of top up-regulated, 10.0% of top up-regulated (padj)\n",paste("Up-regulated :",sum(Normal_Tumor_sort_stage_III$Categories=="Up-regulated"),"Down-regulated :",sum(Normal_Tumor_sort_stage_III$Categories=="Down-regulated"),sep=" "))) 

# Add treshold lines
p2 <- p2 + geom_hline(yintercept=threshold_padj ,linetype = 'dashed') + geom_vline(xintercept=threshold_log2fc_up ,linetype = 'dashed') + geom_vline(xintercept=threshold_log2fc_down ,linetype = 'dashed')

# Selected genes
selected_genes<-rownames(Normal_Tumor_sort_stage_III[Normal_Tumor_sort_stage_III$Categories!="Uncategorized",])
 
# Obtain differential expression numbers
pca_normal_stageIII<-plotPCA(vst_Stage_III_sub[selected_genes,], intgroup="stage_III") + theme_bw() + ggtitle("DE Genes of Stage III") + theme(legend.position='bottom')

Normal_Tumor_sort_stage_III<-Normal_Tumor_sort_stage_III[Normal_Tumor_sort_stage_III$Categories!="Uncategorized",]

# Change histogram plot fill colors by groups
padj_histogram_stage_III<-ggplot(Normal_Tumor_sort_stage_III, aes(x=-log(padj), fill=Categories, color=Categories)) +  geom_histogram(position="identity") + scale_fill_manual(values = c("dodgerblue3", "firebrick3"))  + theme_bw() 
#######################################################################################################################
# FindClusters_resolution
png(filename=paste(output_dir,"Volcano_Plot_Normal_Tumor_Stage3.png",sep=""), width = 28, height = 24, res=600, units = "cm")
	pca_plots<-grid.arrange( p2, padj_histogram_stage_III,pca_normal_stageIII,  ncol = 2)
dev.off()

# Save TSV file with genes from Stage1
write_tsv(Normal_Tumor_sort_stage_III[Normal_Tumor_sort_stage_III$Categories!="Uncategorized",], "/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_III.tsv")
########################################################################################################################













































########################################################################################################################
pca_res <- prcomp(t(assay(vst_stages_sub)))
summ <- summary(pca_res)
summ$importance[2,]


# Compute pca matrix
pca_vst_Stage_I_III    <- data.frame(prcomp(t(assay(vst_stages_sub)))$x)
########################################################################################################################
# First, stage_I_III_III
# Sort table by abs(log2FoldChange) and -log(padj)
df_stage_I_vs_III$up_down<-"uncategorized"
########################################################################################################################
# Select genes and patiens
selected_patients<-colData[colData$stages!="Stage II","patient_id"]
########################################################################################################################
pca_vst_Stage_I_III<-pca_vst_Stage_I_III[selected_patients,]
########################################################################################################################
# Set colnames
pca_vst_Stage_I_III$patient_id<-rownames(pca_vst_Stage_I_III)

# Merge collumns
df_stage_I_vs_III<-merge(pca_vst_Stage_I_III, colData, by = "patient_id")
########################################################################################################################
# Merge collumns
pc1_pc2_Stage_I_vs_III<-ggplot(df_stage_I_vs_III, aes(colour = stages)) + geom_point(aes(x=PC1, y=PC2)) + theme_bw() + theme( legend.position = "none" )  
pc1_pc3_Stage_I_vs_III<-ggplot(df_stage_I_vs_III, aes(colour = stages)) + geom_point(aes(x=PC1, y=PC3)) + theme_bw() + theme( legend.position = "none" )  
pc1_pc4_Stage_I_vs_III<-ggplot(df_stage_I_vs_III, aes(colour = stages)) + geom_point(aes(x=PC1, y=PC4)) + theme_bw() + theme( legend.position = "none" )  
pc2_pc3_Stage_I_vs_III<-ggplot(df_stage_I_vs_III, aes(colour = stages)) + geom_point(aes(x=PC2, y=PC3)) + theme_bw() + theme( legend.position = "none" )  
pc2_pc4_Stage_I_vs_III<-ggplot(df_stage_I_vs_III, aes(colour = stages)) + geom_point(aes(x=PC2, y=PC4)) + theme_bw() + theme( legend.position = "none" )  
pc3_pc4_Stage_I_vs_III<-ggplot(df_stage_I_vs_III, aes(colour = stages)) + geom_point(aes(x=PC3, y=PC4)) + theme_bw() + theme( legend.position = "bottom" )  


# FindClusters_resolution
png(filename=paste(output_dir,"stageI_Components.png",sep=""), width = 18, height = 24, res=600, units = "cm")
	grid.arrange(pc1_pc2_Stage_I_vs_III, pc1_pc3_Stage_I_vs_III, pc1_pc4_Stage_I_vs_III, pc2_pc3_Stage_I_vs_III, pc2_pc4_Stage_I_vs_III, pc3_pc4_Stage_I_vs_III, ncol = 2)  + theme( legend.position = "bottom" )   + theme( legend.position = "none" )  
dev.off()

# FindClusters_resolution
png(filename=paste(output_dir,"pca_stageI_Components.png",sep=""), width = 14, height = 14, res=600, units = "cm")
	scatterplot3d(df_stage_I_vs_III[,c("PC1","PC2","PC3")],color=as.numeric(factor(df_stage_I_vs_III$stages)), pch = 16)
dev.off()

