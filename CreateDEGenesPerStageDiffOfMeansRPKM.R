#######################################################################################################################
# Reload colData from file
# Reload unstranded_data from file
###########################################################################################################################
unstranded_data_file                <- "/home/felipe/Documentos/LungPortal/samples/unstranded_rpkm.tsv"                    #
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
unstranded_data<-unstranded_data[avg_expression_pos$Gene,]

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

# Samples of each stage stored in colData                                                                                             #
sample_stage_I  <-colData[colData$stages=="Stage I","patient_id"]                                                                     #
sample_stage_II <-colData[colData$stages=="Stage II","patient_id"]                                                                    #
sample_stage_III<-colData[colData$stages=="Stage III","patient_id"]                                                                   #

sample_stages_II_III<-colData[colData$stages=="Stage II" | colData$stages=="Stage III","patient_id"]                                  #
sample_stages_I_III<-colData[colData$stages=="Stage I" | colData$stages=="Stage III","patient_id"]                                    #
sample_stages_I_II<-colData[colData$stages=="Stage I" | colData$stages=="Stage II","patient_id"]                                      #
#######################################################################################################################################
df_table_comparisson=rbind(data.frame(Stage_i="sample_stage_I",Stages_ii_and_iii="sample_stages_II_III"),
rbind(data.frame(Stage_i="sample_stage_II",Stages_ii_and_iii="sample_stages_I_III"),
rbind(data.frame(Stage_i="sample_stage_III",Stages_ii_and_iii="sample_stages_I_II"))))
####################################################################################################################
list_of_comparisson=list(sample_stage_I=sample_stage_I,sample_stage_II=sample_stage_II,sample_stage_III=sample_stage_III,sample_stages_II_III=sample_stages_II_III,sample_stages_I_III=sample_stages_I_III,sample_stages_I_II=sample_stages_I_II)
####################################################################################################################
# Create bck for colData_bck
colData_bck<-colData

# for each pair of stage.
for (comparisson_index in rownames(df_table_comparisson))
{	
	# Stages
	Stage_i          <-df_table_comparisson[comparisson_index,"Stage_i"]
	Stages_ii_and_iii<-df_table_comparisson[comparisson_index,"Stages_ii_and_iii"]

	# Take samples of each stage
	Stage_i_samples         =list_of_comparisson[[Stage_i]]
	Stages_ii_and_iii_sample=list_of_comparisson[[Stages_ii_and_iii]]	

	# Take RPKM of genes from samples of each stage
	Stage_i_samples_expr         <-unstranded_data[,Stage_i_samples]
	Stages_ii_and_iii_sample_expr<-unstranded_data[,Stages_ii_and_iii_sample]
	####################################################################################################################
	# folchange=Expr(Stage i)/Expr(Stage ii and II)
	folchange=rowMeans(Stage_i_samples_expr)/rowMeans(Stages_ii_and_iii_sample_expr)

	# log2change
	log2change=log(folchange,2)	

	# log2change data
	log2change_Stage_i=data.frame(gene=names(folchange),log2change=log2change)
	####################################################################################################################	
	# First by padj
	padj_threshold<-1
	log2fc_threshold<-0.58	
	####################################################################################################################
	# First, set category
	# "Unchanged"
	log2change_Stage_i$Category<-"Uncategorized"

	# First, stageI
	log2change_Stage_i[log2change_Stage_i$log2change>=0.58,]<-"Up-regulated"		
  	####################################################################################################################	
	# Save TSV file with genes from Stage3
	write_tsv(log2change_Stage_i, paste(output_dir,"genes_Stage_DiffOfMeansRPKM",Stage_i,".tsv",sep=""))
	####################################################################################################################
	# Selected genes	
	# Obtain differential Category numbers
	selected_genes<-rownames(log2change_Stage_i[which(log2change_Stage_i$Category!="Uncategorized"),])
	 
	# Obtain differential Category numbers
	pca_normal_stages_first_second<-plotPCA(vst_stages_sub[selected_genes,], intgroup=stage_index) + theme_bw() + ggtitle(paste("DE Genes", stage_index))

	# Filter table
	df_stages<-df_stages[which(df_stages$Category!="Uncategorized"),]
	
	# Change histogram plot fill colors by groups
	padj_histogram_stages<-ggplot(df_stages, aes(x=-log(padj), fill=Category, color=Category)) +  geom_histogram(position="identity",bins = 30) + theme_bw() + ggtitle(paste("DE Genes", stage_index))+ theme(legend.position='bottom')
	#######################################################################################################################
	# Remove samples from Stage III and plot again
	pca_normal_stages_first_second<-plotPCA(unstranded_data[selected_genes,], intgroup=stage_index) + theme_bw() + ggtitle(paste("DE Genes ", stage_index,"\n",paste(length(selected_genes), "genes"),sep=""))+ theme(legend.position='bottom')
	#######################################################################################################################
	# FindClusters_resolution
	png(filename=paste(output_dir,"Volcano_Plot_Normal_Tumor_Stage_",stage_index,".png",sep=""), width = 28, height = 24, res=600, units = "cm")
		pca_plots<-grid.arrange( p2, padj_histogram_stages,pca_normal_stages_first_second,  ncol = 2)
	dev.off()		
	#######################################################################################################################	
	# Save TSV file with genes from Stage1
	write_tsv(cbind(data.frame(Gene=rownames(df_stages)),df_stages), paste(output_dir,"DESeq2_selected_genes_Stage_pos_",stage_index,".tsv",sep=""))
}
