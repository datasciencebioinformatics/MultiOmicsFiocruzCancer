#########################################################################################################################################################################
unstranded_rpkm  <-read.table(file = paste("/home/felipe/Documentos/LungPortal/output/","unstranded_rpkm.tsv",sep=""), sep = '\t', header = TRUE,)         #
# Volcano plot of –log10 (p values) vs. log2foldchange. 
# All, log2foldchange tumor,        log2foldchange per stage I, 
#      log2foldchange per stage II, log2foldchange per stage III
# Genes were filtered to include only those that had an RPKM ≥ 1.0, a fold-change of ≥1.5 and a paired t-test value of <0.01. https://pubmed.ncbi.nlm.nih.gov/27884954/
#########################################################################################################################################################################
# Save plots
list_of_plots<-list()
#########################################################################################################################################################################
# First, the log2foldchane tumor/normal samples is used
log2change_tumor_control$Category<-"insignificant"

# First, the log2foldchane tumor/normal samples is used
log2change_tumor_control$Pvalue<-1

# For each genes in the tabe
for (gene in log2change_tumor_control$gene)
{
	# Take p-value
	log2change_tumor_control[gene,"Pvalue"]<-t.test(x=as.numeric(unstranded_rpkm[gene,paired_sample_df$tumor]), y=as.numeric(unstranded_rpkm[gene,paired_sample_df$normal]), paired = TRUE, alternative = "two.sided")$p.value	
}
# FRD 
log2change_tumor_control$FDR<-p.adjust(log2change_tumor_control$Pvalue, method="BH")

# Categorize genes if log2foldchange >= threshold_tumor
log2change_tumor_control[intersect(which(log2change_tumor_control$FDR<=0.05), which(log2change_tumor_control$log2change>=threshold_tumor)),"Category"]<-paste("Tumor genes", sep="")
#########################################################################################################################################################################
# Create volcano plot
p1 <- ggplot(log2change_tumor_control, aes(log2change, -log(FDR,10))) +  geom_point(size = 2/5) +  theme_bw()
	
# Adding color to differentially expressed genes (DEGs)
p1 <- ggplot(log2change_tumor_control, aes(log2change, -log(FDR),color = Category)) + geom_point(size = 2/5,aes(color = Category))  +
  xlab("log2FoldChange") + 
  ylab("-log10(padj)") +
  scale_color_manual(values = c("black", "red")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_bw() + ggtitle(paste("Paired t-test, RPKM of paired tumor/normal samples\nlog2foldchange >=",threshold_tumor, " and FRD <= 0.05", sep="")) + guides(fill="none")

list_of_plots[[tumor]]<-p1
#########################################################################################################################################################################
# for each pair of stage.
for (comparisson_index in rownames(df_table_comparisson))
{	
	# Stages
	Stage_i          <-df_table_comparisson[comparisson_index,"Stage_i"]
	Stages_ii_and_iii<-df_table_comparisson[comparisson_index,"Stages_ii_and_iii"]	
	
	# Take gens of corresponding stage
	DE_genes        <- list_selected_genes[[Stage_i]]
	
	# Take samples of each stage
	Stage_i_samples         =list_of_comparisson[[Stage_i]]
	Stages_ii_and_iii_sample=list_of_comparisson[[Stages_ii_and_iii]]	
	
	# Take RPKM of genes from samples of each stage
	Stage_i_samples_expr         <-na.omit(unstranded_rpkm[,Stage_i_samples])
	Stages_ii_and_iii_sample_expr<-na.omit(unstranded_rpkm[,Stages_ii_and_iii_sample])
	####################################################################################################################
	# Log2foldchange
	LOG_CONSTANT=0.001
	log2change=rowMeans(log(Stage_i_samples_expr+LOG_CONSTANT,2))-rowMeans(log(Stages_ii_and_iii_sample_expr+LOG_CONSTANT,2))
	
	# log2change data
	log2change_tumor_control=na.omit(data.frame(gene=names(log2change),log2change=log2change))
	####################################################################################################################		
	# First, the log2foldchane tumor/normal samples is used
	log2change_tumor_control$Category<-"insignificant"
	
	# First, the log2foldchane tumor/normal samples is used
	log2change_tumor_control$Pvalue<-1
	
	# For each genes in the tabe
	for (gene in log2change_tumor_control$gene)
	{
		# Take p-value
		log2change_tumor_control[gene,"Pvalue"]<-t.test(x=as.numeric(unstranded_rpkm[gene,Stage_i_samples]), y=as.numeric(unstranded_rpkm[gene,Stages_ii_and_iii_sample]), paired = FALSE, alternative = "two.sided")$p.value	
	}
	# FRD 
	log2change_tumor_control$FDR<-p.adjust(log2change_tumor_control$Pvalue, method="BH")

	# Categorize genes if log2foldchange >= threshold_tumor
	log2change_tumor_control[DE_genes,"Category"]<-paste("Tumor genes per stage", sep="")

	# Create volcano plot
	pn <- ggplot(log2change_tumor_control, aes(log2change, -log(FDR,10))) +  geom_point(size = 2/5) +  theme_bw()
		
	# Adding color to differentially expressed genes (DEGs)
	pn <- ggplot(log2change_tumor_control, aes(log2change, -log(FDR),color = Category)) + geom_point(size = 2/5,aes(color = Category))  +
	  xlab("log2FoldChange") + 
	  ylab("-log10(padj)") +
	  scale_color_manual(values = c("black", "red")) +
	  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_bw() + ggtitle(paste("Unpaired t-test, RPKM in ",gsub("sample", "", gsub("_", " ", Stage_i)), " vs. ohter two stages",sep="")) + guides(fill="none")		

	list_of_plots[[Stage_i]]<-pn
}

# Second, the log2foldchane tumor/normal per stage I
# Third, the log2foldchane tumor/normal per stage II
# Fourth, the log2foldchane tumor/normal per stage II
# Save the log2foldchange per stage in a list (from Pipeline_SquamousCellCarcinoma_CreateDEGenesPerStageMeansFromPairedV2.R)
