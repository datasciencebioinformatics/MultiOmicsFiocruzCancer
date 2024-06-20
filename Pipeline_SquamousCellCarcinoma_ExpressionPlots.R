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
sample_stage_I_tumor  <-samples_tumor[samples_tumor$stages=="Stage I",c("tissue_type","patient_id","stages")]                                #
sample_stage_I_normal  <-samples_normal[samples_normal$stages=="Stage I",c("tissue_type","patient_id","stages")]                             #
sample_stage_II_tumor  <-samples_tumor[samples_tumor$stages=="Stage II",c("tissue_type","patient_id","stages")]                              #
sample_stage_II_normal  <-samples_normal[samples_normal$stages=="Stage II",c("tissue_type","patient_id","stages")]                           #
sample_stage_III_tumor  <-samples_tumor[samples_tumor$stages=="Stage III",c("tissue_type","patient_id","stages")]                            #
sample_stage_III_normal  <-samples_normal[samples_normal$stages=="Stage III",c("tissue_type","patient_id","stages")]                         #

# Select samples
sample_stage_all_samples<-rbind(sample_stage_I_tumor,sample_stage_I_normal,sample_stage_II_tumor,sample_stage_II_normal,sample_stage_III_tumor,sample_stage_III_normal)

# Samples
unstranded_data_samples<-melt(t(unstranded_data_filter[,sample_stage_all_samples$patient_id]))

# colnames(unstranded_data_samples)
colnames(unstranded_data_samples)<-c("patient_id","gene_id","RPKM")

# unstranded_data_samples with sample_stage_all_samples
unstranded_data_samples<-merge(unstranded_data_samples,sample_stage_all_samples,by="patient_id")

# colnames
colnames(unstranded_data_samples)<-c("patient_id","gene_id","RPKM","tissue_type","stages","tissue_type","stages_x")

# Change box plot line colors by groups
p<-ggplot(unstranded_data_samples, aes(x=dose, y=len, color=dose)) +   geom_boxplot()
