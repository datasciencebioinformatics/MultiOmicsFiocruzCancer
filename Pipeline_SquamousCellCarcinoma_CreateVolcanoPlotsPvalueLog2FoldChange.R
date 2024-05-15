# Volcano plot of –log10 (p values) vs. log2foldchange. 
# All, log2foldchange tumor,        log2foldchange per stage I, 
#      log2foldchange per stage II, log2foldchange per stage III
# Genes were filtered to include only those that had an RPKM ≥ 1.0, a fold-change of ≥1.5 and a paired t-test value of <0.01. https://pubmed.ncbi.nlm.nih.gov/27884954/
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
	
# The thresholds
threshold_padj<-min(-log(log2change_tumor_control[log2change_tumor_control$Category!="Uncategorized","FDR"]))
threshold_log2fc_up<-min(log2change_tumor_control[log2change_tumor_control$Category=="Up-regulated","log2change"])
threshold_log2fc_down<-max(log2change_tumor_control[log2change_tumor_control$Category=="Down-regulated","log2change"])

# Adding color to differentially expressed genes (DEGs)
p1 <- ggplot(log2change_tumor_control, aes(log2change, -log(FDR),color = Category)) + geom_point(size = 2/5,aes(color = Category))  +
  xlab("log2FoldChange") + 
  ylab("-log10(padj)") +
  scale_color_manual(values = c("black", "red")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_bw() + ggtitle(paste("Paired t-test, RPKM of paired tumor/normal samples\nlog2foldchange >=",threshold_tumor, " and FRD <= 0.05", sep="")) + guides(fill="none")
#########################################################################################################################################################################


# Second, the log2foldchane tumor/normal per stage I
# Third, the log2foldchane tumor/normal per stage II
# Fourth, the log2foldchane tumor/normal per stage II
# Save the log2foldchange per stage in a list (from Pipeline_SquamousCellCarcinoma_CreateDEGenesPerStageMeansFromPairedV2.R)
