########################################################################################################################################
interactome_data_stage_I<-unique(interactome_data_stage_I[,c("Gene1","Gene2")])
interactome_data_stage_II<-unique(interactome_data_stage_II[,c("Gene1","Gene2")])
interactome_data_stage_III<-unique(interactome_data_stage_III[,c("Gene1","Gene2")])
########################################################################################################################################
stage_I_genes_factor  <-factor(c(interactome_data_stage_I$Gene1,interactome_data_stage_I$Gene2),level=unique(c(interactome_data_stage_I$Gene1,interactome_data_stage_I$Gene2)))
stage_II_genes_factor <-factor(c(interactome_data_stage_II$Gene1,interactome_data_stage_II$Gene2),level=unique(c(interactome_data_stage_II$Gene1,interactome_data_stage_II$Gene2)))
stage_III_genes_factor<-factor(c(interactome_data_stage_III$Gene1,interactome_data_stage_III$Gene2),level=unique(c(interactome_data_stage_III$Gene1,interactome_data_stage_III$Gene2)))
   
df_stageI_connectivity   <-unique(data.frame(Conectivity=table(stage_I_genes_factor)))
df_stageII_connectivity  <-unique(data.frame(Conectivity=table(stage_II_genes_factor)))
df_stageIII_connectivity <-unique(data.frame(Conectivity=table(stage_III_genes_factor)))
########################################################################################################################################
colnames(df_stageI_connectivity)<-c("Gene","Conectivity")
colnames(df_stageII_connectivity)<-c("Gene","Conectivity")
colnames(df_stageIII_connectivity)<-c("Gene","Conectivity")
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
########################################################################################################################################
# Save TSV file with genes from Stage1
write_tsv(df_stageI_connectivity, paste(output_dir,"df_stageI_connectivity_I",".tsv",sep=""))
write_tsv(df_stageII_connectivity, paste(output_dir,"df_stageII_connectivity_II",".tsv",sep=""))
write_tsv(df_stageIII_connectivity, paste(output_dir,"df_stageIII_connectivity_II",".tsv",sep=""))

# Save TSV file with genes from Stage1
write_tsv(interactome_data_stage_I, paste(output_dir,"df_stageI_interactome",".tsv",sep=""))
write_tsv(interactome_data_stage_II, paste(output_dir,"df_stageII_interactome",".tsv",sep=""))
write_tsv(interactome_data_stage_III, paste(output_dir,"df_stageIII_interactome",".tsv",sep=""))
########################################################################################################################################
cat(print(paste("\nNº of vertex/Nº/Entropy of edges, co-expression network for Stage I: ",  paste(length(unique(c(interactome_data_stage_I$Gene1,interactome_data_stage_I$Gene2))),dim(unique(interactome_data_stage_I))[1],round(Entropy_stage_I_value_Carels,4),sep="/"),sep="")),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
cat(print(paste("\nNº of vertex/Nº/Entropy of edges, co-expression network for Stage II: ", paste(length(unique(c(interactome_data_stage_II$Gene1,interactome_data_stage_II$Gene2))),dim(unique(interactome_data_stage_II))[1],round(Entropy_stage_II_value_Carels,4),sep="/"),sep="")),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
cat(print(paste("\nNº of vertex/Nº/Entropy of edges, co-expression network for Stage III: ",paste(length(unique(c(interactome_data_stage_III$Gene1,interactome_data_stage_III$Gene2))),dim(unique(interactome_data_stage_III))[1],round(Entropy_stage_III_value_Carels,4),sep="/"),sep="")),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
