#######################################################################################################################
# Reload colData from file
# Reload unstranded_data from file
###########################################################################################################################
unstranded_data_file                <- "/home/felipe/Documentos/LungPortal/samples/unstranded.tsv"                        #
merged_data_patient_info_file       <- "/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv"                  #
colData_file                        <- "/home/felipe/Documentos/LungPortal/samples/colData.tsv"                           #
avg_expression_pos_file             <- "/home/felipe/Documentos/LungPortal/samples/log2fc_expression_pos.tsv"             #
###########################################################################################################################
unstranded_data                    <-read.table(file = unstranded_data_file, sep = '\t', header = TRUE,fill=TRUE)         #
merged_data_patient_info_data      <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)#
colData_data                       <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)                 #
avg_expression_pos                 <-read.table(file = avg_expression_pos_file, sep = '\t', header = TRUE,fill=TRUE)      #
rownames(colData)                  <-colData$patient_id                                                                   #
###########################################################################################################################
# Filter to only positive tumor/normal samples.
#unstranded_data<-unstranded_data[avg_expression_pos$Gene,]

#omit NA values from vector
unstranded_data <- na.omit(unstranded_data)
########################################################################################################################
# A panel to analyse differential Category comparing samples of each stage against all others stages.
########################################################################################################################
# Set colData
colData$stage_I   <- "Stages_II_III"
colData$stage_II  <- "Stages_I_III"
colData$stage_III <- "Stages_I_II"

# Each stage
colData$stage_I[which(colData$stages=="Stage I")]<-"Stage I"
colData$stage_II[which(colData$stages=="Stage II")]<-"Stage II"
colData$stage_III[which(colData$stages=="Stage III")]<-"Stage III"

# Vector with each stage
stages_str<-c("stage_I","stage_II","stage_III")
####################################################################################################################
# Create bck for colData_bck
colData_bck<-colData

# for each pair of stage.
for (stage_index in stages_str)
{	
	# Run DESeq2
	dds_stages <- DESeqDataSetFromMatrix(countData = unstranded_data[,colData$patient_id], colData=colData, design = as.formula(paste("~ age_at_index +  gender + ",stage_index))  )

	# Run DESeq2
	dds_stages <- DESeq(dds_stages)

	print(resultsNames(dds_stages)[4])		

	# Df s6tages I
	df_stages<-data.frame(results(dds_stages,name=resultsNames(dds_stages)[4]))
	####################################################################################################################
	# Run varianceStabilizingTransformation
	vst_stages_sub<-varianceStabilizingTransformation(dds_stages, blind = TRUE, fitType = "parametric")
	####################################################################################################################
	# Remove NA rows
	df_stages<-na.omit(df_stages)	
	####################################################################################################################
	# First by padj
	padj_threshold<-1
	log2fc_threshold<-0.58
	
	# First, stageI
	#df_stages<-df_stages[df_stages$padj<=padj_threshold,]
	####################################################################################################################
	# First, set category
	# "Unchanged"
	df_stages$Category<-"Uncategorized"
	
	# First, stageI
#	df_stages[intersect(which(df_stages$log2FoldChange >= 0),which(abs(df_stages$log2FoldChange)>=log2fc_threshold)),"Category"] <-"Up-regulated"
#	df_stages[intersect(which(df_stages$log2FoldChange <0  ),which(abs(df_stages$log2FoldChange)>=log2fc_threshold)),"Category"] <-"Down-regulated"
	df_stages[intersect(which(df_stages$log2FoldChange >= 0),which(df_stages$log2FoldChange>=log2fc_threshold)),"Category"] <-"Up-regulated"

	# df_stages
	df_stages<-df_stages[rownames(df_stages) %in% avg_expression_pos$Gene,]
  	####################################################################################################################	
	# Save TSV file with genes from Stage3
	write_tsv(cbind(data.frame(Gene=rownames(df_stages)),df_stages), paste(output_dir,"genes_Stage_",stage_index,".tsv",sep=""))
	####################################################################################################################
	# Create volcano plot
	p1 <- ggplot(df_stages, aes(log2FoldChange, -log(padj,2))) +  geom_point(size = 2/5) +  theme_bw()
	
	# The thresholds
	threshold_padj<-min(-log(df_stages[df_stages$Category!="Uncategorized","padj"]))
	threshold_log2fc_up<-min(df_stages[df_stages$Category=="Up-regulated","log2FoldChange"])
	threshold_log2fc_down<-max(df_stages[df_stages$Category=="Down-regulated","log2FoldChange"])

	# Adding color to differentially expressed genes (DEGs)
	p2 <- ggplot(df_stages, aes(log2FoldChange, -log(padj),color = Category)) + geom_point(size = 2/5,aes(color = Category))  +
	  xlab("log2FoldChange") + 
	  ylab("-log(padj)") +
	  scale_color_manual(values = c("gray50", "firebrick3")) +
	  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_bw() + ggtitle(paste("DE Genes ", stage_index,"\n",resultsNames(dds_stages)[4],sep=""))
	
	# Add treshold lines
	p2 <- p2 + geom_hline(yintercept=threshold_padj ,linetype = 'dashed') + geom_vline(xintercept=threshold_log2fc_up ,linetype = 'dashed') + geom_vline(xintercept=threshold_log2fc_down ,linetype = 'dashed')
	
	# Selected genes	
	# Obtain differential Category numbers
	selected_genes<-rownames(df_stages[which(df_stages$Category!="Uncategorized"),])
	 
	# Obtain differential Category numbers
	pca_normal_stages_first_second<-plotPCA(vst_stages_sub[selected_genes,], intgroup=stage_index) + theme_bw() + ggtitle(paste("DE Genes", stage_index))

	# Filter table
	df_stages<-df_stages[which(df_stages$Category!="Uncategorized"),]
	
	# Change histogram plot fill colors by groups
	padj_histogram_stages<-ggplot(df_stages, aes(x=-log(padj), fill=Category, color=Category)) +  geom_histogram(position="identity",bins = 30) + theme_bw() + ggtitle(paste("DE Genes", stage_index))+ theme(legend.position='bottom')
	#######################################################################################################################
	# Remove samples from Stage III and plot again
	pca_normal_stages_first_second<-plotPCA(vst_stages_sub[selected_genes,], intgroup=stage_index) + theme_bw() + ggtitle(paste("DE Genes ", stage_index,"\n",paste(length(selected_genes), "genes"),sep=""))+ theme(legend.position='bottom')
	#######################################################################################################################
	# FindClusters_resolution
	png(filename=paste(output_dir,"Volcano_Plot_Normal_Tumor_Stage_",stage_index,".png",sep=""), width = 28, height = 24, res=600, units = "cm")
		pca_plots<-grid.arrange( p2, padj_histogram_stages,pca_normal_stages_first_second,  ncol = 2)
	dev.off()		
	#######################################################################################################################	
	# Save TSV file with genes from Stage1
	write_tsv(cbind(data.frame(Gene=rownames(df_stages)),df_stages), paste(output_dir,"selected_genes_Stage_pos_",stage_index,".tsv",sep=""))
}
