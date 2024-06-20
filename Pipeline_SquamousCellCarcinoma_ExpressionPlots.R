# Take genes per each stage
genes_stage_I<-gsub(" ", "", unique(unique(strsplit(x=paste(unique(df_count_terms_selected$Genes_Stage_I),collapse=","),split=",",fixed=T))[[1]]))
genes_stage_II<-gsub(" ", "", unique(unique(strsplit(x=paste(unique(df_count_terms_selected$Genes_Stage_II),collapse=","),split=",",fixed=T))[[1]]))
genes_stage_III<-gsub(" ", "", unique(unique(strsplit(x=paste(unique(df_count_terms_selected$Genes_Stage_III),collapse=","),split=",",fixed=T))[[1]]))

##################################################################################################################################################
# Remove empty line
genes_stage_I<-genes_stage_I[genes_stage_I!=""]
genes_stage_II<-genes_stage_II[genes_stage_II!=""]
genes_stage_III<-genes_stage_III[genes_stage_III!=""]

# Select only tumor metadata
colDta_tumor<-colData[colData$tissue_type=="Tumor",]

# Select only tumor samples id's
tumor_samples<-colDta_tumor$patient_id

samples_normal<-colData[paired_sample_df$normal,]
samples_tumor<-colData[paired_sample_df$tumor,]

# Select samples of each stage stored in colData                                                                                             #
sample_stage_I  <-colDta_tumor[samples_normal$stages=="Stage I",]                                                                  #
sample_stage_II <-colDta_tumor[colDta_tumor$stages=="Stage II",]                                                                 #
sample_stage_III <-colDta_tumor[colDta_tumor$stages=="Stage III",]                                                                 #
################################################################################################################################################
