#######################################################################################################################################
g<-graph_from_data_frame(interactome_data, directed = TRUE, vertices = NULL)
gI<-graph_from_data_frame(interactome_dataI, directed = TRUE, vertices = NULL)
gII<-graph_from_data_frame(interactome_dataII, directed = TRUE, vertices = NULL)

interactions_stage_1<-plot(g, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "Network for genes of stage I")
#######################################################################################################################################
# Number of genes
g <- barabasi.game(n=length(genes$gene),m=dim(interactome_data)[1]/length(genes$gene), directed = FALSE, start.graph =g)
as_data_frame(g, what = c("edges"))
interactions_random_1<-plot(g, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "Random network for stage I")
#######################################################################################################################################
N=ceiling((length(genes$gene)+length(genesI$gene)+length(genesII$gene))/3)
#######################################################################################################################################
# entropy_bootstrapping_stage_values
entropy_bootstrapping_random_stage_all<-c()

# Repeat 1000 times
for (bootstrapping in 1:1000)
{
	print(bootstrapping)
	print(Sys.time()) # get start time
	########################################################################################################################################
	# Generate random graph
	g_stage_all_random <- barabasi.game(n=N,m=5, directed = FALSE)
	########################################################################################################################################
	g_stage_all_random<-as_data_frame(g_stage_all_random, what = c("edges"))
	########################################################################################################################################	
	df_random_connectivity_stage_all   <-unique(data.frame(Conectivity=table(c(g_stage_all_random$from,g_stage_all_random$to))))	
	########################################################################################################################################
	colnames(df_random_connectivity_stage_all)<-c("Gene","Conectivity")	
	########################################################################################################################################
	# Table for the calculation of entropy
	df_entropy_calulation_random_stage_all   <-data.frame(table(df_random_connectivity_stage_all$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
	
	# Rename colnames
	colnames(df_entropy_calulation_random_stage_all)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
		
	# Calculate p(k)
	df_entropy_calulation_random_stage_all$p_k<-df_entropy_calulation_random_stage_all$count/sum(df_entropy_calulation_random_stage_all$count)	
	
	# Calculate log2(p(k))
	df_entropy_calulation_random_stage_all$log2_pk<-log(df_entropy_calulation_random_stage_all$p_k,2)
		
	# Calculate p(k)*log2(p(k))
	df_entropy_calulation_random_stage_all$p_k_mult_log2_pk<-df_entropy_calulation_random_stage_all$p_k*df_entropy_calulation_random_stage_all$log2_pk
	
	# Caclulate entropy value
	Entropy_stage_random_value_stage_all    <-abs(sum(df_entropy_calulation_random_stage_all$p_k_mult_log2_pk))

	# Random entropy 
	entropy_bootstrapping_random_stage_all<-c(entropy_bootstrapping_random_stage_all,Entropy_stage_random_value_stage_all)		
}
# Save stages
df_enropy_stage_all  <-data.frame(1:1000,entropy=entropy_bootstrapping_random_stage_all,Method="Conforte")

# Create plot
plot_enropy_stage_all  <-ggplot(df_enropy_stage_all, aes(x=entropy))  + geom_histogram() +  geom_segment(aes(x=Entropy_stage_I_value_Carels, y=200, xend=Entropy_stage_I_value_Carels, yend=0), arrow = arrow(length=unit(0.5, 'cm'))) +  geom_segment(aes(x=Entropy_stage_II_value_Carels, y=200, xend=Entropy_stage_II_value_Carels, yend=0), arrow = arrow(length=unit(0.5, 'cm'))) +  geom_segment(aes(x=Entropy_stage_III_value_Carels, y=225, xend=Entropy_stage_III_value_Carels, yend=0), arrow = arrow(length=unit(0.5, 'cm')))  + theme_bw() + annotate(geom="text", x=Entropy_stage_I_value_Carels, y=210, label="Stage I")+ annotate(geom="text", x=Entropy_stage_II_value_Carels, y=210, label="Stage II")+ annotate(geom="text", x=Entropy_stage_III_value_Carels, y=230, label="Stage III")

# FindClusters_resolution
png(filename=paste(output_dir,"Entropy_","all_.png",sep=""), width = 16, height = 16, res=600, units = "cm")
	plot_enropy_stage_all + ggtitle("Entropy bootstrapping (1000x random gene sets)")
dev.off()
