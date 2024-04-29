####################################################################################################################
unstranded_data_file                <- "/home/felipe/Documentos/LungPortal/samples/unstranded.tsv"                 #
unstranded_data                     <-read.table(file = unstranded_data_file, sep = '\t', header = TRUE,fill=TRUE) #
####################################################################################################################interactome_data_stage_I
merge_interactome_data<-rbind(data.frame(Gene1=interactome_data_stage_I_clean$Gene1,Gene2=interactome_data_stage_I_clean$Gene2,Stage="Stage I"),
data.frame(Gene1=interactome_data_stage_II_clean$Gene1,Gene2=interactome_data_stage_II_clean$Gene2,Stage="Stage II"),
data.frame(Gene1=interactome_data_stage_III_clean$Gene1,Gene2=interactome_data_stage_III_clean$Gene2,Stage="Stage III"))

# merge_interactome_data
merge_interactome_data<-unique(merge_interactome_data[,c("Gene1","Gene2")])

# Number of genes
n_of_interactions<-ceiling((dim(interactome_data_stage_I)[1]+ dim(interactome_data_stage_II)[1]+dim(interactome_data_stage_III)[1])/3)

# entropy_bootstrapping_stage_values
entropy_bootstrapping_Carels<-c()
entropy_bootstrapping_Shannon<-c()

# Repeat 1000 times
for (bootstrapping in 1:1000)
{
  	# Random genes 
	random_interactions_Stage_all<-sample( 1:dim(interactome_data)[1], n_of_interactions, replace = TRUE, prob = NULL)  

	# Take random interactions
	random_interactions<-interactome_data[random_interactions_Stage_all,]	

	# Colnames
	colnames(random_interactions)<-c("Gene1","Gene2")
	
	# PPI counts
	interactome_data_stage_all   <-unique(data.frame(Conectivity=table(c(random_interactions$Gene1,random_interactions$Gene2))))

	# Rename names
	colnames(interactome_data_stage_all)<-c("Gene","Conectivity")

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
	
	# Caclulate entropy value
	Entropy_value_Shannon_stage_random<-Entropy(df_stageI_connectivity$Conectivity, base=exp(2))

	entropy_bootstrapping_Carels<-c(entropy_bootstrapping_Carels,Entropy_stage_random_value_Carels)
	entropy_bootstrapping_Shannon<-c(entropy_bootstrapping_Shannon,Entropy_value_Shannon_stage_random)	
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
####################################################################################################################
round(sum(entropy_stage_I<=entropy_bootstrapping_stage_values)/1000,3)
round(sum(entropy_stage_II<=entropy_bootstrapping_stage_values)/1000,3)
round(sum(entropy_stage_III<=entropy_bootstrapping_stage_values)/1000,3)
####################################################################################################################

