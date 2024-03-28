####################################################################################################################
# Filter down colData_bck:
# First,  only stages I and III
# Second, only stages II and III
# Third,  only stages II and III
####################################################################################################################
# Aalysed pairs
df_stage_pairs<-data.frame(First=c("Stage I","Stage I","Stage II"),Second=c("Stage II","Stage III","Stage III"))
####################################################################################################################
# Create bck for colData_bck
colData_bck<-colData

# for each pair of stage.
for (stage_pair in rownames(df_stage_pairs))
{
	# First and second stages
	first_stage<-df_stage_pairs[stage_pair,"First"]
	second_stage<-df_stage_pairs[stage_pair,"Second"]
	
	# Create bck for colData_bck
	colData_bck<-colData

	# Filter col data
	colData_bck<-colData_bck[union(which(colData_bck$stages==first_stage), which(colData_bck$stages==second_stage)),]

	# Re-factor stages
	colData_bck$stages<-factor(colData_bck$stages)
			
	# Run DESeq2
	dds_stages <- DESeqDataSetFromMatrix(countData = unstranded_data[,colData_bck$patient_id], colData=colData_bck, design = ~ age_at_index +  gender +stages  )

	# Run DESeq2
	dds_stages <- DESeq(dds_stages)

	print(resultsNames(dds_stages)[4])		

	# Df s6tages I
	df_stages<-data.frame(results(dds_stages,name=resultsNames(dds_stages)[4]))

	####################################################################################################################
	# Run varianceStabilizingTransformation
	vst_stages_sub<-varianceStabilizingTransformation(dds_stages, blind = TRUE, fitType = "parametric")
	####################################################################################################################
	# First, stage O
	# First, stageI
	# Sort table by abs(log2FoldChange) and -log(padj)
	down_regulated_df_stage<-df_stages[which(df_stages$log2FoldChange<=0), ]
	up_regulated_df_stages<-df_stages[which(df_stages$log2FoldChange>0), ]
	####################################################################################################################
	# Sort table by abs(log2FoldChange) and -log(padj)
	Normal_Tumor_up_sort_df_stages<- up_regulated_df_stages[order(up_regulated_df_stages$padj), ]
	Normal_Tumor_down_sort_df_stages<- down_regulated_df_stage[order(down_regulated_df_stage$padj), ]
	####################################################################################################################
	# Remove NA rows
	# First, stageI
	Normal_Tumor_up_sort_df_stages<-na.omit(Normal_Tumor_up_sort_df_stages)
	Normal_Tumor_down_sort_df_stages<-na.omit(Normal_Tumor_down_sort_df_stages)
	####################################################################################################################
	# Field for top 10 percent of sorted sample
	# First, stageI
	Normal_Tumor_up_sort_df_stages$Normal_Tumor_sort_percent<-FALSE
	Normal_Tumor_down_sort_df_stages$Normal_Tumor_sort_percent<-FALSE
	Normal_Tumor_up_sort_df_stages$up_down<-"up"
	Normal_Tumor_down_sort_df_stages$up_down<-"down"
	####################################################################################################################
	# First, stageI
	# Field for top 10 percent of sorted sample
	Normal_Tumor_up_sort_df_stages[1:(dim(Normal_Tumor_up_sort_df_stages)[1]*0.10),"Normal_Tumor_sort_percent"]<-TRUE
	Normal_Tumor_down_sort_df_stages[1:(dim(Normal_Tumor_down_sort_df_stages)[1]*0.10),"Normal_Tumor_sort_percent"]<-TRUE
	####################################################################################################################
	Normal_Tumor_sort_stages<-rbind(Normal_Tumor_up_sort_df_stages,Normal_Tumor_down_sort_df_stages)
	####################################################################################################################
	# First, stageI
	# "Unchanged"
	Normal_Tumor_sort_stages$Expression<-0
	####################################################################################################################
	# Set expression up
	# First, stageI
	Normal_Tumor_sort_stages[intersect(which(Normal_Tumor_sort_stages$Normal_Tumor_sort_percent), which(Normal_Tumor_sort_stages$log2FoldChange < 0)),"Expression"]<--1
	Normal_Tumor_sort_stages[intersect(which(Normal_Tumor_sort_stages$Normal_Tumor_sort_percent), which(Normal_Tumor_sort_stages$log2FoldChange >= 0)),"Expression"]<-1
	####################################################################################################################
	# First, stageI
	Normal_Tumor_sort_stages$Categories<-""
	Normal_Tumor_sort_stages[which(Normal_Tumor_sort_stages$Expression==0),"Categories"]<-"Uncategorized"
	Normal_Tumor_sort_stages[which(Normal_Tumor_sort_stages$Expression==1),"Categories"]<-"Up-regulated"
	Normal_Tumor_sort_stages[which(Normal_Tumor_sort_stages$Expression==-1),"Categories"]<-"Down-regulated"
	####################################################################################################################
	# Save TSV file with genes from Stage3
	write_tsv(Normal_Tumor_sort_stages, paste(output_dir,"genes_Stage",stage_pair,".tsv",sep=""))
	####################################################################################################################
	# Create volcano plot
	p1 <- ggplot(Normal_Tumor_sort_stages, aes(log2FoldChange, -log(padj,2))) +  geom_point(size = 2/5) +  theme_bw()
	
	# The thresholds
	threshold_padj<-min(-log(Normal_Tumor_sort_stages[Normal_Tumor_sort_stages$Categories!="Uncategorized","padj"]))
	threshold_log2fc_up<-min(Normal_Tumor_sort_stages[Normal_Tumor_sort_stages$Categories=="Up-regulated","log2FoldChange"])
	#threshold_log2fc_down<-max(Normal_Tumor_sort_stages[Normal_Tumor_sort_stages$Categories=="Down-regulated","log2FoldChange"])

	# Adding color to differentially expressed genes (DEGs)
	p2 <- ggplot(Normal_Tumor_sort_stages, aes(log2FoldChange, -log(padj),color = Categories)) + geom_point(size = 2/5,aes(color = Categories))  +
	  xlab(expression("log2FoldChange")) + 
	  ylab(expression("-log(padj)")) +
	  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
	  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_bw() + ggtitle(paste("DE Genes", first_stage, "vs.",second_stage))
	
	# Add treshold lines
	p2 <- p2 + geom_hline(yintercept=threshold_padj ,linetype = 'dashed') + geom_vline(xintercept=threshold_log2fc_up ,linetype = 'dashed') + geom_vline(xintercept=threshold_log2fc_down ,linetype = 'dashed')
	
	# Selected genes
	
	# Obtain differential expression numbers
	selected_genes<-rownames(Normal_Tumor_sort_stages[which(Normal_Tumor_sort_stages[which(Normal_Tumor_sort_stages$Categories=="Up-regulated"),"Normal_Tumor_sort_percent"]),])		
	 
	# Obtain differential expression numbers
	pca_normal_stages_first_second<-plotPCA(vst_stages_sub[selected_genes,], intgroup="stages") + theme_bw() + ggtitle(paste("DE Genes", first_stage, "vs.",second_stage))
	
	Normal_Tumor_sort_stages<-Normal_Tumor_sort_stages[selected_genes,]
	
	# Change histogram plot fill colors by groups
	padj_histogram_stages<-ggplot(Normal_Tumor_sort_stages, aes(x=-log(padj), fill=Categories, color=Categories)) +  geom_histogram(position="identity") + theme_bw()  + theme_bw() + ggtitle(paste("DE Genes", first_stage, " vs.",second_stage,"\n",paste(sum(Normal_Tumor_sort_stages$Categories=="Up-regulated"), "genes\n10% of genes sorted by padj\nselected up-regulated genes"),sep=""))+ theme(legend.position='bottom')
	#######################################################################################################################
	# Remove samples from Stage III and plot again
	pca_normal_stages_first_second<-plotPCA(vst_stages_sub[selected_genes,], intgroup="stages") + theme_bw() + ggtitle(paste("DE Genes", first_stage, " vs.",second_stage,"\n",paste(sum(Normal_Tumor_sort_stages$Categories=="Up-regulated"), "genes\n20% of genes sorted by padj\nselected up-regulated genes"),sep=""))+ theme(legend.position='bottom')
	#######################################################################################################################
	# FindClusters_resolution
	png(filename=paste(output_dir,"Volcano_Plot_Normal_Tumor_Stage",stage_pair,".png",sep=""), width = 28, height = 24, res=600, units = "cm")
		pca_plots<-grid.arrange( p2, padj_histogram_stages,pca_normal_stages_first_second,  ncol = 2)
	dev.off()		
	
	# Save TSV file with genes from Stage1
	# Save TSV file with genes from Stage1
	Normal_Tumor_sort_stages$Gene<-rownames(Normal_Tumor_sort_stages)
	write_tsv(Normal_Tumor_sort_stages[Normal_Tumor_sort_stages$Categories=="Up-regulated",], paste(output_dir,"selected_genes_Stage_pos",stage_pair,".tsv",sep=""))
}
