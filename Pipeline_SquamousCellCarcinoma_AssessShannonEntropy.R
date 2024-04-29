#######################################################################################################################################
# Load interactome data
interactome_data <-read.table(file = interactome_file, sep = '\t', header = FALSE,fill=TRUE)         

# Rename collumns 
colnames(interactome_data)<-c("Gene1","Gene2")
#######################################################################################################################################
# Filter tables to keep only the gene entries that are listed in the EnsemblToUniprotKBconversionList
interactome_data<-interactome_data[interactome_data$Gene1 %in% EnsemblToUniprotKBconversionList_data$SYMBOL,]
interactome_data<-interactome_data[interactome_data$Gene2 %in% EnsemblToUniprotKBconversionList_data$SYMBOL,]

# Create a table for id conversion gene_id and gene_symbol for the genes in the interactome data
gene1_conversion<-merge(interactome_data,EnsemblToUniprotKBconversionList_data,by.x="Gene1", by.y="SYMBOL",all.x=TRUE,all.y=FALSE)
gene_conversion<-merge(gene1_conversion,EnsemblToUniprotKBconversionList_data,by.x="Gene2", by.y="SYMBOL",all.x=TRUE,all.y=FALSE)

# Keep only the collumns of interest- 
# interactome_data : interactome with converted ids
interactome_data<-unique(gene_conversion[,3:4])

# Rename interactome_data collumns
colnames(interactome_data)<-c("Gene1","Gene2")

# Filter tables to keep only the gene entries that are listed in the EnsemblToUniprotKBconversionList
interactome_data<-interactome_data[interactome_data$Gene1 %in% genes_ids,]
interactome_data<-interactome_data[interactome_data$Gene2 %in% genes_ids,]
#######################################################################################################################################
# merge_interactome_data
merge_interactome_data<-unique(merge_interactome_data[,c("Gene1","Gene2")])

# Number of genes
n_of_interactions<-ceiling((length(genes_interactome_stage_I)+ length(genes_interactome_stage_II)+length(genes_interactome_stage_II))/3)

# entropy_bootstrapping_stage_values
entropy_bootstrapping_Carels<-c()

# Repeat 1000 times
for (bootstrapping in 1:1000)
{
	# Store genes stage I, II and III
	# Vectors to store gene ids from each stage
	genes_id_vector_random<-sample( genes_ids, n_of_interactions, replace = TRUE, prob = NULL)  

	# Colnames
	colnames(random_interactions)<-c("Gene1","Gene2")

	# If at least one of the genes in the pair are in the interactome
	interactome_data_stage_random<-rbind(interactome_data[interactome_data$Gene1 %in% genes_id_vector_random,],
	interactome_data[interactome_data$Gene2 %in% genes_id_vector_random,])

	# Clean the tables
	for (gene_pair_index in rownames(interactome_data_stage_random))
	{
	    # interactome_data_stage
	    pair_gene_id_I <-interactome_data_stage_random[gene_pair_index,"Gene1"]
	    pair_gene_id_II<-interactome_data_stage_random[gene_pair_index,"Gene2"]
	
	    # If both genes are in the list of genes_ids
	    if( (pair_gene_id_I %in% genes_ids) &&  (pair_gene_id_II %in% genes_ids) )
	    {
	      # Re-order gene ids
	      if(pair_gene_id_II<pair_gene_id_I)
	      {
	        interactome_data_stage_random[gene_pair_index,"Gene1"]<-pair_gene_id_II
	        interactome_data_stage_random[gene_pair_index,"Gene2"]<-pair_gene_id_I     
	      }
	      # Re-order gene ids
	      if(pair_gene_id_II==pair_gene_id_I)
	      {
	        interactome_data_stage_random[gene_pair_index,"Gene2"]<-"REPEAT"
	      }
	    }
	}
	# Take unique values
	interactome_data_stage_random<-unique(interactome_data_stage_random)

	# df_stageI_connectivity
	df_random_connectivity   <-unique(data.frame(Conectivity=table(c(interactome_data_stage_random$Gene1,interactome_data_stage_random$Gene2))))
		
	# PPI counts
	interactome_data_stage_all   <-unique(data.frame(Conectivity=table(c(interactome_data_stage_random$Gene1,interactome_data_stage_random$Gene2))))

	# Rename names
	colnames(interactome_data_stage_all)<-c("Gene","Conectivity")

	# Remove REPEAT
	interactome_data_stage_all<-interactome_data_stage_all[interactome_data_stage_all$Gene!="REPEAT",]

	# Table for the calculation of entropy
	df_entropy_calulation_random   <-data.frame(table(interactome_data_stage_all$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
	
	# Rename colnames
	colnames(df_entropy_calulation_random)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
	
	# Calculate p(k)
	df_entropy_calulation_random$p_k<-df_entropy_calulation_random$count/sum(df_entropy_calulation_random$count)
	
	# Calculate log2(p(k))
	df_entropy_calulation_random$log2_pk<-log(df_entropy_calulation_random$p_k,2)
	
	# Calculate p(k)*log2(p(k))
	df_entropy_calulation_random$p_k_mult_log2_pk<-df_entropy_calulation_random$p_k*df_entropy_calulation_random$log2_pk
	
	# Caclulate entropy value
	Entropy_stage_random_value_Carels  <-abs(sum(df_entropy_calulation_random$p_k_mult_log2_pk))

	# Store in the vector
	entropy_bootstrapping_Carels<-c(entropy_bootstrapping_Carels,Entropy_stage_random_value_Carels)
}
# Save stages
df_enropy_stage_all  <-data.frame(1:1000,entropy=entropy_bootstrapping_Carels,Method="Conforte")

# Create plot
plot_enropy_stage_all<-ggplot(df_enropy_stage_all, aes(x=entropy))  + geom_histogram() 

# Histogram overlaid with kernel density curve
plot_enropy_stage_all     <-plot_enropy_stage_all +  geom_segment(aes(x=Entropy_stage_I_value_Carels, y=200, xend=Entropy_stage_I_value_Carels, yend=0), arrow = arrow(length=unit(0.5, 'cm'))) +   geom_segment(aes(x=Entropy_stage_II_value_Carels, y=200, xend=Entropy_stage_II_value_Carels, yend=0), arrow = arrow(length=unit(0.5, 'cm'))) +  geom_segment(aes(x=Entropy_stage_III_value_Carels, y=200, xend=Entropy_stage_III_value_Carels, yend=0), arrow = arrow(length=unit(0.5, 'cm'))) + theme_bw()
plot_enropy_stage_all     <-plot_enropy_stage_all + annotate("text", x = Entropy_stage_I_value_Carels, y = 205, label = "Stage I")
plot_enropy_stage_all     <-plot_enropy_stage_all + annotate("text", x = Entropy_stage_II_value_Carels, y = 205, label = "Stage II") 
plot_enropy_stage_all     <- plot_enropy_stage_all + annotate("text", x = Entropy_stage_III_value_Carels, y = 205, label = "Stage III")

# FindClusters_resolution
png(filename=paste(output_dir,"Entropy_","all_.png",sep=""), width = 16, height = 16, res=600, units = "cm")
	plot_enropy_stage_all + ggtitle("Entropy bootstrapping (1000x random gene sets)")
dev.off()

sum(Entropy_stage_I_value_Carels<=entropy_bootstrapping_Carels)/1000
sum(Entropy_stage_II_value_Carels<=entropy_bootstrapping_Carels)/1000
sum(Entropy_stage_III_value_Carels<=entropy_bootstrapping_Carels)/1000
