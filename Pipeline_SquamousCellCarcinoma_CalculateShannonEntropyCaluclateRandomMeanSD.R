# Gene table
genes_Stage_I       <-read.table(file = file_genes_Stage_I, sep = '\t', header = TRUE,fill=TRUE)         
genes_Stage_II      <-read.table(file = file_genes_Stage_II, sep = '\t', header = TRUE,fill=TRUE)
genes_Stage_III     <-read.table(file = file_genes_Stage_III, sep = '\t', header = TRUE,fill=TRUE)   

# Stage-specific stage I
# Stage-specific stage II
# Stage-specific stage III
# Save the values of entropy for stage I
entropy_stage_I<-c()
entropy_stage_II<-c()
entropy_stage_III<-c()

for (random in 1:1000)
{
  random_genes_I   <-sample(genes_ids, dim(genes_Stage_I)[1])
  random_genes_II  <-sample(genes_ids, dim(genes_Stage_II)[1])
  random_genes_III <-sample(genes_ids, dim(genes_Stage_III)[1])

  full_interactome_stage_I<- data.frame(expand.grid.unique(x = genes_interactome_stage_I, y = genes_interactome_stage_I,include.equals=FALSE))
  full_interactome_stage_II<- data.frame(expand.grid.unique(x = genes_interactome_stage_II, y = genes_interactome_stage_II,include.equals=FALSE))
  full_interactome_stage_III<- data.frame(expand.grid.unique(x = genes_interactome_stage_III, y = genes_interactome_stage_III,include.equals=FALSE))

    # set colnames
  colnames(full_interactome_stage_I)<-c("Gene1","Gene2")
  colnames(full_interactome_stage_II)<-c("Gene1","Gene2")
  colnames(full_interactome_stage_III)<-c("Gene1","Gene2")
  
  rownames(full_interactome_stage_I)<-paste(full_interactome_stage_I$Gene1,full_interactome_stage_I$Gene2,sep="-")
  rownames(full_interactome_stage_II)<-paste(full_interactome_stage_II$Gene1,full_interactome_stage_II$Gene2,sep="-")
  rownames(full_interactome_stage_III)<-paste(full_interactome_stage_III$Gene1,full_interactome_stage_III$Gene2,sep="-")
  #######################################################################################################
  interactome_stage_I  <-interactome_data_inv[which(rownames(interactome_data_inv) %in% rownames(full_interactome_stage_I)),]
  interactome_stage_II <-interactome_data_inv[which(rownames(interactome_data_inv) %in% rownames(full_interactome_stage_II)),]
  interactome_stage_III<-interactome_data_inv[which(rownames(interactome_data_inv) %in% rownames(full_interactome_stage_III)),]
  ########################################################################################################################################
  stage_I_genes_factor  <-factor(c(interactome_stage_I$Gene1,interactome_stage_I$Gene2),level=unique(c(interactome_stage_I$Gene1,interactome_stage_I$Gene2)))
  stage_II_genes_factor <-factor(c(interactome_stage_II$Gene1,interactome_stage_II$Gene2),level=unique(c(interactome_stage_II$Gene1,interactome_stage_II$Gene2)))
  stage_III_genes_factor<-factor(c(interactome_stage_III$Gene1,interactome_stage_III$Gene2),level=unique(c(interactome_stage_III$Gene1,interactome_stage_III$Gene2)))
     
  df_stageI_connectivity   <-unique(data.frame(Conectivity=table(stage_I_genes_factor)))
  df_stageII_connectivity  <-unique(data.frame(Conectivity=table(stage_II_genes_factor)))
  df_stageIII_connectivity <-unique(data.frame(Conectivity=table(stage_III_genes_factor)))
  
  #df_stageI_connectivity   <-unique(data.frame(Conectivity=table(c(interactome_stage_I$Gene1,interactome_stage_I$Gene2))))
  #df_stageII_connectivity  <-unique(data.frame(Conectivity=table(c(interactome_stage_II$Gene1,interactome_stage_II$Gene2))))
  #df_stageIII_connectivity <-unique(data.frame(Conectivity=table(c(interactome_stage_III$Gene1,interactome_stage_III$Gene2))))
  ########################################################################################################################################
  # If dim length equal to zero
  if(dim(df_stageI_connectivity)[1]==0)
  {
      df_stageI_connectivity<-data.frame(Conectivity.stage_I_genes_factor=c("REMOVE"),Conectivity.Freq=0)
  }
  # If dim length equal to zero
  if(dim(df_stageII_connectivity)[1]==0)
  {
      df_stageII_connectivity<-data.frame(Conectivity.stage_II_genes_factor=c("REMOVE"),Conectivity.Freq=0)
  }  
  # If dim length equal to zero
  if(dim(df_stageIII_connectivity)[1]==0)
  {
      df_stageIII_connectivity<-data.frame(Conectivity.stage_III_genes_factor=c("REMOVE"),Conectivity.Freq=0)
  }    
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

  entropy_stage_I<-c(entropy_stage_I,Entropy_stage_I_value_Carels)
  entropy_stage_II<-c(entropy_stage_II,Entropy_stage_II_value_Carels)
  entropy_stage_III<-c(entropy_stage_III,Entropy_stage_III_value_Carels)
}
  






