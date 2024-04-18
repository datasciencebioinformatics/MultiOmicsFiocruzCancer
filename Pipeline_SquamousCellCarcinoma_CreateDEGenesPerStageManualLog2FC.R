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
# Writing mtcars data                                                                                                     #
df_log2foldchange_StageI_file = "/home/felipe/Documentos/LungPortal/samples/df_log2foldchange_stageI.tsv"                 #
df_log2foldchange_StageII_file = "/home/felipe/Documentos/LungPortal/samples/df_log2foldchange_stageII.tsv"               #
df_log2foldchange_StageIII_file = "/home/felipe/Documentos/LungPortal/samples/df_log2foldchange_stageIII.tsv"             #
                                                                                                                          #
df_log2foldchange_StageI_data   <-read.table(file = df_log2foldchange_StageI_file, sep = '\t', header = TRUE,fill=TRUE)   #
df_log2foldchange_StageII_data  <-read.table(file = df_log2foldchange_StageII_file, sep = '\t', header = TRUE,fill=TRUE)  #
df_log2foldchange_StageIII_data <-read.table(file = df_log2foldchange_StageIII_file, sep = '\t', header = TRUE,fill=TRUE) #
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

# Save data.frames
list_genes_per_stage<-list(stage_I=df_log2foldchange_StageI_data,stage_II=df_log2foldchange_StageII_data,stage_III=df_log2foldchange_StageIII_data)
####################################################################################################################
# Create bck for colData_bck
colData_bck<-colData

# for each pair of stage.
for (stage_index in stages_str)
{		
	####################################################################################################################
	# Df s6tages I
	df_stages<-list_genes_per_stage[[stage_index]]
	####################################################################################################################
  	colnames(df_stages)<-c("Gene","log2foldchange")	
	####################################################################################################################	
	# First by padj
	padj_threshold<-1
	log2fc_threshold<-0
	
	# First, stage
	df_stages<-df_stages[df_stages$log2foldchange>=log2fc_threshold,]
	####################################################################################################################
	# First, set category
	# "Unchanged"
	df_stages$Category<-"Uncategorized"
	
	# First, stageI
	df_stages[which(df_stages$log2foldchange>=log2fc_threshold),"Category"] <-"Up-regulated"

	# df_stages
	#df_stages<-df_stages[rownames(df_stages) %in% avg_expression_pos$Gene,]
  ####################################################################################################################	
	# Save TSV file with genes from Stage3
	write_tsv(cbind(data.frame(Gene=rownames(df_stages)),df_stages), paste(output_dir,"genes_Stage_manual_log2folchange",stage_index,".tsv",sep=""))
	####################################################################################################################
	# Selected genes	
	# Obtain differential Category numbers
	selected_genes<-rownames(df_stages[which(df_stages$Category!="Uncategorized"),])
  ####################################################################################################################
	# Run DESeq2
	dds_stages <- DESeqDataSetFromMatrix(countData = na.omit(unstranded_data[selected_genes,colData$patient_id]) , colData=colData, design = as.formula(paste("~ age_at_index +  gender + ",stage_index))  )

	# Run DESeq2
	dds_stages <- DESeq(dds_stages)

	print(resultsNames(dds_stages)[4])		

	# Df s6tages I
	df_results<-data.frame(results(dds_stages,name=resultsNames(dds_stages)[4]))
	
	# Run varianceStabilizingTransformation
	vst_stages_sub<-varianceStabilizingTransformation(dds_stages, blind = TRUE, fitType = "parametric")
	 
	# Obtain differential Category numbers
	pca_normal_stages_first_second<-plotPCA(vst_stages_sub, intgroup=stage_index) + theme_bw() + ggtitle(paste("DE Genes", stage_index))

	# Filter table
	df_stages<-df_stages[which(df_stages$Category!="Uncategorized"),]	
	#######################################################################################################################
	# Remove samples from Stage III and plot again
	pca_normal_stages_first_second<-plotPCA(vst_stages_sub, intgroup=stage_index) + theme_bw() + ggtitle(paste("DE Genes ", stage_index,"\n",paste(length(selected_genes), "genes"),sep=""))+ theme(legend.position='bottom')
	#######################################################################################################################
	# FindClusters_resolution
	png(filename=paste(output_dir,"PCA_Plot_Normal_Tumor_Stage_",stage_index,".png",sep=""), width = 16, height = 16, res=600, units = "cm")
		plotPCA(vst_stages_sub, intgroup=stage_index) + theme_bw() + ggtitle(paste("DE Genes ", stage_index,"\n",paste(length(selected_genes), "genes"),sep=""))+ theme(legend.position='bottom')
	dev.off()		
	#######################################################################################################################	
	# Save TSV file with genes from Stage1
	write_tsv(cbind(data.frame(Gene=rownames(df_stages)),df_stages), paste(output_dir,"Manual_selected_genes_Stage_pos_",stage_index,".tsv",sep=""))
}
