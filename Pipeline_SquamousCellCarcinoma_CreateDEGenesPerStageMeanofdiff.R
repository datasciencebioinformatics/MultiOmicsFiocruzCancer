# Number of genes per stage
genes_id_file_stageI  <- "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageI/genes_id.txt"
genes_id_file_stageII <- "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageII/genes_id.txt"
genes_id_file_stageIII<- "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageIII/genes_id.txt"

# Read table
genes_id_data_stageI   <-read.table(file = genes_id_file_stageI, sep = '\t', header = FALSE,fill=TRUE)  
genes_id_data_stageII  <-read.table(file = genes_id_file_stageII, sep = '\t', header = FALSE,fill=TRUE)  
genes_id_data_stageIII <-read.table(file = genes_id_file_stageIII, sep = '\t', header = FALSE,fill=TRUE)  
#####################################################################################################################
# Reading the contents of TSV file using read_tsv() method
cases_files_stageI<- "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageI/cases_files.txt"

# Load metadata table
df_cases_files_stageI     <-read.table(file = cases_files_stageI, sep = '\t', header = FALSE,fill=TRUE)  

# Data frame of results
df_results_stageI<-data.frame(Gene=genes_id_data_stageI$V1)

# vector for the ids
cases_ids_1<-c()

# For each df_cases_files_stageI
for (file in df_cases_files_stageI[,1])
{
	print(file)
	df_cases_filesI_stageI     <-read.table(file = file, sep = '\t', header = TRUE,fill=TRUE)  
					
	# Re-order
	df_cases_filesI_stageI<-data.frame(Gene=genes_id_data_stageI$V1,Diff=df_cases_filesI_stageI[genes_id_data_stageI$V1,])
	
	# Merge all tables
	df_results_stageI<-merge(df_results_stageI,df_cases_filesI_stageI,by="Gene",all = TRUE)
	
	# Save cases_ids
	cases_ids_1<-c(cases_ids_1,colnames(df_results_stageI))		
}
#####################################################################################################################
# Reading the contents of TSV file using read_tsv() method
cases_files_stageII<- "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageII/cases_files.txt"

# Load metadata table
df_cases_files_stageII     <-read.table(file = cases_files_stageII, sep = '\t', header = FALSE,fill=TRUE)  

# Data frame of results
df_results_stageII<-data.frame(Gene=genes_id_data_stageII$V1)

# vector for the ids
cases_ids_2<-c()

# For each df_cases_files_stageI
for (file in df_cases_files_stageII[,1])
{
	print(file)
	df_cases_filesII_stageII     <-read.table(file = file, sep = '\t', header = TRUE,fill=TRUE)  
					
	# Re-order
	df_cases_filesII_stageII<-data.frame(Gene=genes_id_data_stageII$V1,Diff=df_cases_filesII_stageII[genes_id_data_stageII$V1,])
	
	# Merge all tables
	df_results_stageII<-merge(df_results_stageII,df_cases_filesII_stageII,by="Gene",all = TRUE)
	
	# Save cases_ids
	cases_ids_2<-c(cases_ids_2,colnames(df_results_stageII))	
}
#####################################################################################################################
# Reading the contents of TSV file using read_tsv() method
cases_files_stageIII<- "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageIII/cases_files.txt"

# Load metadata table
df_cases_files_stageIII     <-read.table(file = cases_files_stageIII, sep = '\t', header = FALSE,fill=TRUE)  

# Data frame of results
df_results_stageIII<-data.frame(Gene=genes_id_data_stageIII$V1)

# vector for the ids
cases_ids_3<-c()

# For each df_cases_files_stageI
for (file in df_cases_files_stageIII[,1])
{
	print(file)
	df_cases_filesIII_stageIII     <-read.table(file = file, sep = '\t', header = TRUE,fill=TRUE)  
					
	# Re-order
	df_cases_filesIII_stageIII<-data.frame(Gene=genes_id_data_stageIII$V1,Diff=df_cases_filesIII_stageIII[genes_id_data_stageIII$V1,])
	
	# Merge all tables
	df_results_stageIII<-merge(df_results_stageIII,df_cases_filesIII_stageIII,by="Gene",all = TRUE)
	
	# Save cases_ids
	cases_ids_3<-c(cases_ids_3,colnames(df_cases_filesIII_stageIII))
}
#####################################################################################################################
rownames(df_results_stageI)<-df_results_stageI$Gene
rownames(df_results_stageII)<-df_results_stageII$Gene
rownames(df_results_stageIII)<-df_results_stageIII$Gene

df_results_stageI <-df_results_stageI[,-1]
df_results_stageII<-df_results_stageII[,-1]
df_results_stageIII<-df_results_stageIII[,-1]

df_results_stageI <-data.frame(Avg_log2folchange=rowMeans(df_results_stageI,na.rm=TRUE),gene=rownames(df_results_stageI))
df_results_stageII<-data.frame(Avg_log2folchange=rowMeans(df_results_stageII,na.rm=TRUE),gene=rownames(df_results_stageII))
df_results_stageIII<-data.frame(Avg_log2folchange=rowMeans(df_results_stageIII,na.rm=TRUE),gene=rownames(df_results_stageIII))
#####################################################################################################################
genes_results_stageI  <- merge_interactome_gene_symbol[which(merge_interactome_gene_symbol$gene_name %in% df_results_stageI$gene),"gene_name"]
genes_results_stageII <- merge_interactome_gene_symbol[which(merge_interactome_gene_symbol$gene_name %in% df_results_stageII$gene),"gene_name"]
genes_results_stageIII<-  merge_interactome_gene_symbol[which(merge_interactome_gene_symbol$gene_name %in% df_results_stageIII$gene),"gene_name"]
#####################################################################################################################
df_results_stageI<-df_results_stageI[genes_results_stageI,]
df_results_stageII<-df_results_stageII[genes_results_stageII,]
df_results_stageIII<-df_results_stageIII[genes_results_stageIII,]
#####################################################################################################################
# Save data.frames
df_results_stageI$gene_id  <-merge_interactome_gene_symbol[which(merge_interactome_gene_symbol$gene_name %in% df_results_stageI$gene),"gene_id"]
df_results_stageII$gene_id  <-merge_interactome_gene_symbol[which(merge_interactome_gene_symbol$gene_name %in% df_results_stageII$gene),"gene_id"]
df_results_stageIII$gene_id  <-merge_interactome_gene_symbol[which(merge_interactome_gene_symbol$gene_name %in% df_results_stageIII$gene),"gene_id"]
#####################################################################################################################
write_tsv(df_results_stageI, paste(output_dir,"Meanofdiff_selected_genes_Stage_pos_","I",".tsv",sep=""))
write_tsv(df_results_stageII, paste(output_dir,"Meanofdiff_selected_genes_Stage_pos_","II",".tsv",sep=""))
write_tsv(df_results_stageIII, paste(output_dir,"Meanofdiff_selected_genes_Stage_pos_","III",".tsv",sep=""))
#####################################################################################################################
