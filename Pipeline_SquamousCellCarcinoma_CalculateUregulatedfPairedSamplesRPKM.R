###########################################################################################################################
# Data.frame to store results for control samples
df_diff_expression<-data.frame(Gene=c(),log2foldchange=c(),stage=c())

# set rownames
rownames(df_diff_expression)<-df_diff_expression$Gene

# For each stages
for (stage in c("Stage I","Stage II","Stage III"))
{
	# Subset annotation per stage
	colData_stage_i<-colData[colData$stages==stage,]
		
	# take normal and tumor	
	normal_sample<-paired_sample_df[,"normal"]
	tumor_sample <-paired_sample_df[,"tumor"]    
		
	# Take patients from stage i, normal and tumor	
	patient_stage_i_normal<-colData_stage_i[colData_stage_i$patient_id %in% normal_sample,"patient_id"]
	patient_stage_i_tumor<-colData_stage_i[colData_stage_i$patient_id %in% tumor_sample,"patient_id"]
	
	# Store results for the case
	case_results_normal<-data.frame(case=norm_counts[,patient_stage_i_normal])
	case_results_tumor <-data.frame(case=norm_counts[,patient_stage_i_tumor])	
	
	# Calculate differential expression
	diff_expression<-log(rowMeans(case_results_tumor),2)-log(rowMeans(case_results_normal),2)
	
	# Store results
	df_stage_i_resuts<-data.frame(Gene=names(diff_expression),log2foldchange=diff_expression,stage=stage)

	# Concatenate 	
	df_diff_expression<-rbind(df_diff_expression,df_stage_i_resuts)
}
# remove omit
df_diff_expression<-na.omit(df_diff_expression)

# remove infinite
df_diff_expression <- df_diff_expression[is.finite(df_diff_expression$log2foldchange),]
###########################################################################################################################
sum(df_diff_expression[df_diff_expression$stage=="Stage I","log2foldchange"]>=3.00)
sum(df_diff_expression[df_diff_expression$stage=="Stage II","log2foldchange"]>=3.00)
sum(df_diff_expression[df_diff_expression$stage=="Stage III","log2foldchange"]>=3.00)

# df stages
df_stage_I   <-df_diff_expression[which(df_diff_expression[df_diff_expression$stage=="Stage I","log2foldchange"]>=3.00),]
df_stage_II  <-df_diff_expression[which(df_diff_expression[df_diff_expression$stage=="Stage II","log2foldchange"]>=3.00),]
df_stage_III <-df_diff_expression[which(df_diff_expression[df_diff_expression$stage=="Stage III","log2foldchange"]>=3.00),]

# Write table
write.table(df_stage_I, file = "/home/felipe/Documentos/LungPortal/samples/log2foldchange_stageI_pos", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(df_stage_II, file = "/home/felipe/Documentos/LungPortal/samples/log2foldchange_stageII_pos", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(df_stage_III, file = "/home/felipe/Documentos/LungPortal/samples/log2foldchange_stageIII_pos", sep = "\t", row.names = TRUE, col.names = TRUE)
###########################################################################################################################
library(ggfortify) 
library(ggplot2)

# For each stages
for (stage in c("Stage I","Stage II","Stage III"))
{  
	# Subset annotation per stage
	colData_stage_i<-colData[colData$stages==stage,]
	
	# take normal and tumor	
	normal_sample<-paired_sample_df[,"normal"]
	tumor_sample <-paired_sample_df[,"tumor"]    
	
	# Take patients from stage i, normal and tumor	
	patient_stage_i_normal<-colData_stage_i[colData_stage_i$patient_id %in% normal_sample,"patient_id"]
	patient_stage_i_tumor<-colData_stage_i[colData_stage_i$patient_id %in% tumor_sample,"patient_id"]
	
	# Selected genes	
	# Obtain differential Category numbers
	selected_genes<-df_stage_I$Gene
	
	pca_res <- prcomp(t(unstranded_data[selected_genes,]), scale. = TRUE)
	dt_pca <- data.frame('Stages' = colData[colnames(unstranded_data[selected_genes,]),"stages"], pca_res$x[,1:2])		
	
	
	# FindClusters_resolution
	png(filename=paste(output_dir,"PCA_tumor_control_",stage,".png",sep=""), width = 16, height = 16, res=600, units = "cm")
	print(print(ggplot2::autoplot(pca_res, data=colData[colnames(unstranded_data[selected_genes,]),], colour="stages", frame=FALSE, frame.type="t") + xlim(-0.15,0.15) + ylim(-0.15,0.15) + theme_bw() + ggtitle(paste("DE Genes ", Stage_i,"\n",paste(length(selected_genes), "genes"),paste("\n",length(patient_stage_i_tumor)," tumor/control paired samples",sep=""),sep=""))+ theme(legend.position='bottom')))
	dev.off()
}
##########################################################################################################################
