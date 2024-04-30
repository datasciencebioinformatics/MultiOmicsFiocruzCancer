#######################################################################################################################################
# Number of genes
n_of_interactions<-ceiling((length(genes_interactome_stage_I)+ length(genes_interactome_stage_II)+length(genes_interactome_stage_II))/3)

# entropy_bootstrapping_stage_values
entropy_bootstrapping_Stage_I_Carels<-c()
entropy_bootstrapping_Stage_II_Carels<-c()
entropy_bootstrapping_Stage_III_Carels<-c()

# Repeat 1000 times
for (bootstrapping in 1:1000)
{
	print(bootstrapping)
	
	# Store genes stage I, II and III
	# Vectors to store gene ids from each stage
	genes_id_vector_random<-sample( genes_ids, n_of_interactions, replace = TRUE, prob = NULL)  

	# Vectors to store gene ids from each stage
	df_correlation_net_stage_I<-t(data.frame(na.omit(unstranded_data_filter[genes_id_vector_random,sample_stage_I])))
	df_correlation_net_stage_II<-t(data.frame(na.omit(unstranded_data_filter[genes_id_vector_random,sample_stage_II])))
	df_correlation_net_stage_III<-t(data.frame(na.omit(unstranded_data_filter[genes_id_vector_random,sample_stage_III])))

	# Merge data
	merge_interactome_data<-rbind(data.frame(Gene1=interactome_data_stage_I$Gene1,Gene2=interactome_data_stage_I$Gene2,Stage="Stage I"),
	data.frame(Gene1=interactome_data_stage_II$Gene1,Gene2=interactome_data_stage_II$Gene2,Stage="Stage II"),
	data.frame(Gene1=interactome_data_stage_III$Gene1,Gene2=interactome_data_stage_III$Gene2,Stage="Stage III"))
	
	# Clean the tables
	for (gene_pair_index in rownames(merge_interactome_data))
	{
	    # interactome_data_stage
	    pair_gene_id_I <-as.vector(merge_interactome_data[gene_pair_index,"Gene1"])
	    pair_gene_id_II<-as.vector(merge_interactome_data[gene_pair_index,"Gene2"])
	
	    # Re-order gene ids
	    if(pair_gene_id_II<pair_gene_id_I)
	    {
	      merge_interactome_data[gene_pair_index,"Gene1"]<-pair_gene_id_II
	      merge_interactome_data[gene_pair_index,"Gene2"]<-pair_gene_id_I     
	    }
	    # Re-order gene ids
	    if(pair_gene_id_II==pair_gene_id_I)
	    {
	      merge_interactome_data[gene_pair_index,"Gene2"]<-"REPEAT"
	    }
	}
	# Take unique values
	merge_interactome_data<-unique(merge_interactome_data)
	########################################################################################################################################
	interactome_data_stage_I<-merge_interactome_data[merge_interactome_data$Stage=="Stage I",]
	interactome_data_stage_II<-merge_interactome_data[merge_interactome_data$Stage=="Stage II",]
	interactome_data_stage_III<-merge_interactome_data[merge_interactome_data$Stage=="Stage III",]
	########################################################################################################################################
	interactome_data_stage_I<-unique(interactome_data_stage_I[,c("Gene1","Gene2")])
	interactome_data_stage_II<-unique(interactome_data_stage_II[,c("Gene1","Gene2")])
	interactome_data_stage_III<-unique(interactome_data_stage_III[,c("Gene1","Gene2")])
	########################################################################################################################################
	df_stageI_connectivity   <-unique(data.frame(Conectivity=table(c(interactome_data_stage_I$Gene1,interactome_data_stage_I$Gene2))))
	df_stageII_connectivity  <-unique(data.frame(Conectivity=table(c(interactome_data_stage_II$Gene1,interactome_data_stage_II$Gene2))))
	df_stageIII_connectivity <-unique(data.frame(Conectivity=table(c(interactome_data_stage_III$Gene1,interactome_data_stage_III$Gene2))))
	########################################################################################################################################
	colnames(df_stageI_connectivity)<-c("Gene","Conectivity")
	colnames(df_stageII_connectivity)<-c("Gene","Conectivity")
	colnames(df_stageIII_connectivity)<-c("Gene","Conectivity")
	########################################################################################################################################
	df_stageI_connectivity<-df_stageI_connectivity[df_stageI_connectivity$Gene!="REPEAT",]
	df_stageII_connectivity<-df_stageII_connectivity[df_stageII_connectivity$Gene!="REPEAT",]
	df_stageIII_connectivity<-df_stageIII_connectivity[df_stageIII_connectivity$Gene!="REPEAT",]
	########################################################################################################################################
	# Table for the calculation of entropy
	df_entropy_calulation_I   <-data.frame(table(df_stageI_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
	df_entropy_calulation_II  <-data.frame(table(df_stageII_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
	df_entropy_calulation_III <-data.frame(table(df_stageIII_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
	
	# Rename colnames
	colnames(df_entropy_calulation_I)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
	colnames(df_entropy_calulation_II)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
	colnames(df_entropy_calulation_III)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
	
	# Calculate p(k)
	df_entropy_calulation_I$p_k<-df_entropy_calulation_I$count/sum(df_entropy_calulation_I$count)
	df_entropy_calulation_II$p_k<-df_entropy_calulation_II$count/sum(df_entropy_calulation_II$count)
	df_entropy_calulation_III$p_k<-df_entropy_calulation_III$count/sum(df_entropy_calulation_III$count)
	
	# Calculate log2(p(k))
	df_entropy_calulation_I$log2_pk<-log(df_entropy_calulation_I$p_k,2)
	df_entropy_calulation_II$log2_pk<-log(df_entropy_calulation_II$p_k,2)
	df_entropy_calulation_III$log2_pk<-log(df_entropy_calulation_III$p_k,2)
	
	# Calculate p(k)*log2(p(k))
	df_entropy_calulation_I$p_k_mult_log2_pk<-df_entropy_calulation_I$p_k*df_entropy_calulation_I$log2_pk
	df_entropy_calulation_II$p_k_mult_log2_pk<-df_entropy_calulation_II$p_k*df_entropy_calulation_II$log2_pk
	df_entropy_calulation_III$p_k_mult_log2_pk<-df_entropy_calulation_III$p_k*df_entropy_calulation_III$log2_pk
	
	# Caclulate entropy value
	Entropy_stage_I_value_Carels  <-abs(sum(df_entropy_calulation_I$p_k_mult_log2_pk))
	Entropy_stage_II_value_Carels <-abs(sum(df_entropy_calulation_II$p_k_mult_log2_pk))
	Entropy_stage_III_value_Carels<-abs(sum(df_entropy_calulation_III$p_k_mult_log2_pk))

	entropy_bootstrapping_Stage_I_Carels<-c(entropy_bootstrapping_Stage_I_Carels,Entropy_stage_I_value_Carels)
	entropy_bootstrapping_Stage_II_Carels<-c(entropy_bootstrapping_Stage_II_Carels,Entropy_stage_II_value_Carels)
	entropy_bootstrapping_Stage_III_Carels<-c(entropy_bootstrapping_Stage_II_Carels,Entropy_stage_III_value_Carels)
		
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
