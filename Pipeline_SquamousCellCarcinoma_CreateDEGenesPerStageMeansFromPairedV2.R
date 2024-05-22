#######################################################################################################################
# Reload colData from file
# Reload unstranded_data from file
###########################################################################################################################
merged_data_patient_info_file       <- "/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv"                  #
colData_file                        <- "/home/felipe/Documentos/LungPortal/samples/colData.tsv"                           #
###########################################################################################################################
unstranded_data                    <-unstranded_data_filter
merged_data_patient_info_data      <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)#
colData_data                       <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)                 #
rownames(colData)                  <-colData$patient_id                                                                   #
###########################################################################################################################
#omit NA values from vector
unstranded_data <- na.omit(unstranded_data)
########################################################################################################################
# A panel to analyse differential Category comparing samples of each stage against all others stages.
########################################################################################################################
# Only tumor samples
colData_tumor<-colData[colData$tissue_type=="Tumor",]

# Set colData
colData_tumor$stage_I   <- "Stages_II_III"
colData_tumor$stage_II  <- "Stages_I_III"
colData_tumor$stage_III <- "Stages_I_II"

# Vector with each stage
stages_str<-c("stage_I","stage_II","stage_III")

# Samples of each stage stored in colData                                                                                             #
sample_stage_I  <-colData_tumor[colData_tumor$stages=="Stage I","patient_id"]                                                                     #
sample_stage_II <-colData_tumor[colData_tumor$stages=="Stage II","patient_id"]                                                                    #
sample_stage_III<-colData_tumor[colData_tumor$stages=="Stage III","patient_id"]                                                                   #

sample_stages_II_III<-colData_tumor[colData_tumor$stages=="Stage II" | colData_tumor$stages=="Stage III","patient_id"]                                  #
sample_stages_I_III<-colData_tumor[colData_tumor$stages=="Stage I" | colData_tumor$stages=="Stage III","patient_id"]                                    #
sample_stages_I_II<-colData_tumor[colData_tumor$stages=="Stage I" | colData_tumor$stages=="Stage II","patient_id"]                                      #
#######################################################################################################################################
df_table_comparisson=rbind(data.frame(Stage_i="sample_stage_I",Stages_ii_and_iii="sample_stages_II_III"),
rbind(data.frame(Stage_i="sample_stage_II",Stages_ii_and_iii="sample_stages_I_III"),
rbind(data.frame(Stage_i="sample_stage_III",Stages_ii_and_iii="sample_stages_I_II"))))
####################################################################################################################
list_of_comparisson=list(sample_stage_I=sample_stage_I,sample_stage_II=sample_stage_II,sample_stage_III=sample_stage_III,sample_stages_II_III=sample_stages_II_III,sample_stages_I_III=sample_stages_I_III,sample_stages_I_II=sample_stages_I_II)
#list_of_genes=list(genes_stage_I=df_stage_I_data$Gene,genes_stage_II=df_stage_II_data$Gene,genes_stage_III=df_stage_III_data$Gene)
list_of_genes=list(genes_stage_I=rownames(na.omit(unstranded_data_filter)),genes_stage_II=rownames(na.omit(unstranded_data_filter)),genes_stage_III=rownames(na.omit(unstranded_data_filter)))
####################################################################################################################
# Create bck for colData_bck
colData_bck<-colData


# for each pair of stage.
for (comparisson_index in rownames(df_table_comparisson))
{	
	# Stages
	Stage_i          <-df_table_comparisson[comparisson_index,"Stage_i"]
	Stages_ii_and_iii<-df_table_comparisson[comparisson_index,"Stages_ii_and_iii"]
	
	print(names(list_of_genes)[as.integer(comparisson_index)])
	
	# Take gens of corresponding stage
	DE_genes        <- list_of_genes[[as.integer(comparisson_index)]]
	
	# Take samples of each stage
	Stage_i_samples         =list_of_comparisson[[Stage_i]]
	Stages_ii_and_iii_sample=list_of_comparisson[[Stages_ii_and_iii]]	
	
	# Take RPKM of genes from samples of each stage
	Stage_i_samples_expr         <-na.omit(unstranded_data[DE_genes,Stage_i_samples])
	Stages_ii_and_iii_sample_expr<-na.omit(unstranded_data[DE_genes,Stages_ii_and_iii_sample])
	####################################################################################################################
	# folchange=Expr(Stage i)/Expr(Stage ii and II)
	#folchange=rowMeans(Stage_i_samples_expr)/rowMeans(Stages_ii_and_iii_sample_expr)	
	#log2change=log(folchange,2)	

	# Log2foldchange
	LOG_CONSTANT=0.001
	log2change=rowMeans(log(Stage_i_samples_expr+LOG_CONSTANT,2))/rowMeans(log(Stages_ii_and_iii_sample_expr+LOG_CONSTANT,2))
	log2change=log( (rowMeans(Stage_i_samples_expr+LOG_CONSTANT)/rowMeans(Stages_ii_and_iii_sample_expr+LOG_CONSTANT)),2)	
	
	# log2change data
	log2change_Stage_i=na.omit(data.frame(gene=names(log2change),log2change=log2change))

	# First, the log2foldchane tumor/normal samples is used
	log2change_Stage_i$Category<-"insignificant"

	# First, the log2foldchane tumor/normal samples is used
	log2change_Stage_i$Pvalue<-1

	# For each genes in the tabe
	for (gene in log2change_Stage_i$gene)
	{
		# Take p-value
		log2change_Stage_i[gene,"Pvalue"]<-t.test(x=as.numeric(unstranded_data[gene,Stage_i_samples]), y=as.numeric(unstranded_data[gene,Stages_ii_and_iii_sample]), paired = FALSE, alternative = "two.sided")$p.value	
	}
	# FRD 
	log2change_Stage_i$FDR<-p.adjust(log2change_Stage_i$Pvalue, method="fdr")

	# Categorize genes if log2foldchange >= threshold_tumor
	log2change_Stage_i[intersect(which(log2change_Stage_i$FDR<=threshold_FDR), which(log2change_Stage_i$log2change>=threshold_stage)),"Category"]<-paste("Per stage genes", sep="")
	#log2change_Stage_i[which(log2change_Stage_i$FDR<=0.05),"Category"]<-paste("Per stage genes", sep="")
	
	# Selected genes based on FDR, log2foldchange
	selected_genes<-log2change_Stage_i[intersect(which(log2change_Stage_i$FDR<=threshold_FDR), which(log2change_Stage_i$log2change>=threshold_stage)),"gene"]	
	####################################################################################################################	
	# Save TSV file with genes from Stage3
	write_tsv(na.omit(log2change_Stage_i[selected_genes,]), paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_Stage_",Stage_i,".tsv",sep=""))			####################################################################################################################		
	cat(print(paste("\nNumber of tumor genes per stage for ",Stage_i, " : ",length(selected_genes))),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)

	#######################################################################################################################################
	write_tsv(log2change_Stage_i, paste(output_dir,gsub('sample_s', 'S' ,Stage_i),"log2change_Stage_i_table.tsv",sep=""))			
	####################################################################################################################	
	# Create volcano plot
	p1 <- ggplot(log2change_Stage_i, aes(log2change, -log(FDR,10))) +  geom_point(size = 2/5) +  theme_bw()
		
	# Adding color to differentially expressed genes (DEGs)
	p1 <- ggplot(log2change_Stage_i, aes(log2change, -log(FDR),color = Category)) + geom_point(size = 2/5,aes(color = Category))  +
	  xlab("log2FoldChange") + 
	  ylab("-log10(padj)") +
	  scale_color_manual(values = c("black", "red")) +
	  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_bw() + ggtitle(paste("Unpaired t-test, RPKM of samples from ", gsub('sample_s', 'S' ,Stage_i), "\nlog2foldchange >=",threshold_tumor, " and FRD <= ",threshold_FDR, sep="")) + guides(fill="none")
	#######################################################################################################################################		
	png(filename=paste(output_dir,"Volcano_plot_",gsub('sample_s', 'S' ,Stage_i),".png",sep=""), width = 16, height = 16, res=600, units = "cm")                                                                                                    #
		p1
	dev.off()
}
