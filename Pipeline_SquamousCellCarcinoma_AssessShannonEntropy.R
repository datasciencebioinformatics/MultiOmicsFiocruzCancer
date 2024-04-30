#######################################################################################################################################
# Number of genes
n_of_interactions<-ceiling((length(genes_interactome_stage_I)+ length(genes_interactome_stage_II)+length(genes_interactome_stage_II))/3)

# entropy_bootstrapping_stage_values
entropy_bootstrapping_random_Carels<-c()

# Repeat 1000 times
for (bootstrapping in 1:1000)
{
	print(bootstrapping)
	print(Sys.time()) # get start time

	# Data frame of df_random_interacions
	df_random_interacions<-data.frame(Gene1=c(),Gene2=c())

	# For each interactions
	for (i in 1:n_of_interactions*(n_of_interactions-1)/2)
	{
		# Store genes stage I, II and III
		# Vectors to store gene ids from each stage
		node_pair_1<-sample( rownames(na.omit(unstranded_data_filter)), 1, replace = TRUE, prob = NULL) 
		node_pair_2<-sample( rownames(na.omit(unstranded_data_filter)), 1, replace = TRUE, prob = NULL) 
		
		# Re-order gene ids
		if(pair_gene_id_II<pair_gene_id_I)
		{
			node_pair_tmp<-node_pair_1
			node_pair_1  <-node_pair_2
			node_pair_2  <-node_pair_1
		}
		# Re-order gene ids
		if(pair_gene_id_II==pair_gene_id_I)
		{
			node_pair_2 <- "REPEAT"
		}	
		# If probability greater 0.7
		if (runif(1, min=0, max=1)>0.1)
		{
			df_random_interacions<-rbind(df_random_interacions,data.frame(Gene1=node_pair_1,Gene2=node_pair_2))
		}	
	}
	colnames(df_random_interacions)<-c("Gene1","Gene2")

	########################################################################################################################################
	df_random_connectivity   <-unique(data.frame(Conectivity=table(c(df_random_interacions$Gene1,df_random_interacions$Gene2))))	
	########################################################################################################################################
	colnames(df_random_connectivity)<-c("Gene","Conectivity")	
	########################################################################################################################################
	# Table for the calculation of entropy
	df_entropy_calulation_random   <-data.frame(table(df_random_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
	
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

	# Random entropy 
	entropy_bootstrapping_random_Carels<-c(entropy_bootstrapping_random_Carels,Entropy_stage_random_value_Carels)		
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
