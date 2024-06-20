# Take genes per each stage
genes_stage_I<-gsub(" ", "", unique(unique(strsplit(x=paste(unique(df_count_terms_selected$Genes_Stage_I),collapse=","),split=",",fixed=T))[[1]]))
genes_stage_II<-gsub(" ", "", unique(unique(strsplit(x=paste(unique(df_count_terms_selected$Genes_Stage_II),collapse=","),split=",",fixed=T))[[1]]))
genes_stage_III<-gsub(" ", "", unique(unique(strsplit(x=paste(unique(df_count_terms_selected$Genes_Stage_III),collapse=","),split=",",fixed=T))[[1]]))
