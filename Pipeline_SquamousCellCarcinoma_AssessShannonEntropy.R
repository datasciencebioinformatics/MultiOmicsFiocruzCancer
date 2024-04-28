####################################################################################################################
unstranded_data_file                <- "/home/felipe/Documentos/LungPortal/samples/unstranded.tsv"                 #
unstranded_data                     <-read.table(file = unstranded_data_file, sep = '\t', header = TRUE,fill=TRUE) #
####################################################################################################################
# Take stageI_list_of_genes
list_of_genes<-unique(merge_interactome_gene_symbol[,c("gene_id","PPI")])

# File path
file_genes_Stage_I   <-   "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_I.tsv"
file_genes_Stage_II   <-  "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_II.tsv"
file_genes_Stage_III   <- "/home/felipe/Documentos/LungPortal/output/DE_GenesPerStageMeansFromPairedUp_unique_stage_III.tsv"

# Gene table
genes_Stage_I       <-read.table(file = file_genes_Stage_I, sep = '\t', header = TRUE,fill=TRUE)         #
genes_Stage_II      <-read.table(file = file_genes_Stage_II, sep = '\t', header = TRUE,fill=TRUE)#
genes_Stage_III     <-read.table(file = file_genes_Stage_III, sep = '\t', header = TRUE,fill=TRUE)

# genes_Stages
genes_Stage_I<-rownames(unstranded_data_filter)[rownames(unstranded_data_filter) %in% genes_Stage_I$gene]
genes_Stage_II<-rownames(unstranded_data_filter)[rownames(unstranded_data_filter) %in% genes_Stage_II$gene]
genes_Stage_III<-rownames(unstranded_data_filter)[rownames(unstranded_data_filter) %in% genes_Stage_III$gene]

# Number of genes
n_of_genes<-ceiling((length(genes_Stage_I)+ length(genes_Stage_II)+length(genes_Stage_III))/3)

# entropy_bootstrapping_stage_values
entropy_bootstrapping_stage_values<-c()

# Repeat 1000 times
for (bootstrapping in 1:1000)
{
  	# Random genes 
	random_genes_Stage_all<-sample(rownames(unstranded_data), n_of_genes, replace = FALSE, prob = NULL)  

	# vector to store all genes 
	genes_id_vector_stage_all<-c()

	# For each gene in stage I
	for (random_gene in random_genes_Stage_all)
	{
	  # Store gene id in the vector
	  genes_id_vector_stage_all<-c(genes_id_vector_stage_all,strsplit(random_gene, split = "\\.")[[1]][1])
	}	

	# If at least one of the genes in the pair are in the interactome
	interactome_data_stage_all<-unique(rbind(interactome_data[interactome_data$Gene1 %in% genes_id_vector_stage_all,],
	interactome_data[interactome_data$Gene2 %in% genes_id_vector_stage_all,]))

	# PPI counts
	interactome_data_stage_all   <-unique(data.frame(Conectivity=table(c(interactome_data_stage_all$Gene1,interactome_data_stage_all$Gene2))))

	# entropy_bootstrapping_samples
	entropy_bootstrapping_stage_values<-c(entropy_bootstrapping_stage_values,round(Entropy(interactome_data_stage_all$Conectivity.Freq, base = 2),3))
}
# Save stages
df_enropy_stage_all  <-data.frame(1:1000,entropy=entropy_bootstrapping_stage_values,stage="Stages")

plot_enropy_stage_all<-ggplot(df_enropy_stage_all, aes(x=entropy))  + geom_histogram() 

# Histogram overlaid with kernel density curve
plot_enropy_stage_all     <-plot_enropy_stage_all +  geom_segment(aes(x=entropy_stage_I, y=200, xend=entropy_stage_I, yend=0), arrow = arrow(length=unit(0.5, 'cm'))) +  geom_segment(aes(x=entropy_stage_II, y=200, xend=entropy_stage_II, yend=0), arrow = arrow(length=unit(0.5, 'cm'))) +  geom_segment(aes(x=entropy_stage_III, y=200, xend=entropy_stage_III, yend=0), arrow = arrow(length=unit(0.5, 'cm'))) 

# FindClusters_resolution
png(filename=paste(output_dir,"Entropy_","all_.png",sep=""), width = 16, height = 16, res=600, units = "cm")
	plot_enropy_stage_all
dev.off()
