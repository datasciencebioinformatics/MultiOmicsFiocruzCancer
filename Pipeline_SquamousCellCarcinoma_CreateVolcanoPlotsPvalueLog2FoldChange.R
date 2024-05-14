# Volcano plot of –log10 (p values) vs. log2foldchange. 
# All, log2foldchange tumor,        log2foldchange per stage I, 
#      log2foldchange per stage II, log2foldchange per stage III

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
