library(readr)
library("xlsx")
library(ggplot2)
library("DESeq2")
##########################################################################################################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output/"                                                             #
#####################################################################################################################
# Reading the contents of TSV file using read_tsv() method
merged_data_patient_info_file<- "/home/felipe/Documentos/LungPortal/samples/merged_data_patient_info.tsv"

# Load metadata table
merged_data_patient_info     <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)   
#####################################################################################################################
# A script to load cancer data base in R
unstranded_file       <- "/home/felipe/Documentos/LungPortal/samples/unstranded.rna_seq.augmented_star_gene_counts.tsv"
colnames_file         <- "/home/felipe/Documentos/LungPortal/samples/header.txt"
rownames_file         <- "/home/felipe/Documentos/LungPortal/samples/gene_ids.txt"

# Load data
unstranded_data<-read.table(file = unstranded_file, sep = '\t', header = FALSE,fill=TRUE)    
colnames_data<-read.table(file = colnames_file, sep = '\t', header = FALSE,fill=TRUE)                                    
rownames_data<-read.table(file = rownames_file, sep = '\t', header = FALSE,fill=TRUE)      

# Set colnames and rownames
rownames(unstranded_data)<-rownames_data[,1]
colnames(unstranded_data)<-colnames_data[,1]

# To check counts and sample_ids
sum(merged_data_patient_info$sample_id %in% colnames(unstranded_data))

# A field to store 
merged_data_patient_info$stages<-merged_data_patient_info$ajcc_pathologic_stage

# Group stages I,II,III and IV
merged_data_patient_info$stages<-gsub("Stage IA", "Stage I", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IB", "Stage I", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIA", "Stage II", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIB", "Stage II", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIIA", "Stage III", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIIB", "Stage III", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IVA", "Stage IV", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IVB", "Stage IV", merged_data_patient_info$stages)

# Filter the table to contain only samples in the dataset
merged_data_patient_info<-merged_data_patient_info[merged_data_patient_info$sample_id %in% colnames(unstranded_data),]

# Save table with samples and metatada
# There are 238 cases, 285 samples, 570 pairs

# Filter collumns that are used for age_at_index, gender, stages, Sample.ID
merged_data_patient_info<-merged_data_patient_info[,c("case_id","sample_id","age_at_index","gender","Sample.Type","stages")]

# Rename collumns
colnames(merged_data_patient_info)<-c("case_id","sample_id","age_at_index","gender","tumor_normal","stages")

#####################################################################################################################
# Set colData
colData<-unique(merged_data_patient_info[,c("sample_id","age_at_index","gender","tumor_normal","stages")])

# Set colnames
rownames(colData)<-colData$sample_id
#####################################################################################################################
# To DO:
# Filter stages with N>30
merged_data_patient_info<-merged_data_patient_info[merged_data_patient_info$stage %in%  names(which(table(merged_data_patient_info$stage)>30)),]

# Filter DESeq2 steps
#####################################################################################################################
# Organize how to send to Carles
write_tsv(unstranded_data, "/home/felipe/Documentos/LungPortal/samples/unstranded.rna_seq.gene_counts.tsv")
write_tsv(colData, "/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv")
#####################################################################################################################
# Add columns for patient ID
colData$patient_id<-paste("patient_",1:length(colData$sample_id),sep="")

# Set colnames
rownames(colData)<-colData$sample_id

# Set rownames
colnames(unstranded_data) <- colData[colnames(unstranded_data),"patient_id"]

# Set colnames
rownames(colData)<-colData$patient_id
#####################################################################################################################
# Organize how to send to Carles
write_tsv(unstranded_data, "/home/felipe/Documentos/LungPortal/samples/unstranded.rna_seq.gene_counts.tsv")
write_tsv(colData, "/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv")
#####################################################################################################################
# Set collumn to factor
#####################################################################################################################
# Set collumn to factor
colData$sample_id<-factor(colData$sample_id)
colData$tumor_normal<-factor(colData$tumor_normal)
colData$stages<-factor(colData$stages)
colData$gender<-factor(colData$gender)
colData$age_range<-factor(cut(colData$age_at_index, breaks = c(0, 25, 50,75,100 )))

# Creat DEseq element from unstranded_data and merged_data_patient_info
dds <- DESeqDataSetFromMatrix(countData = unstranded_data, colData=colData[colnames(unstranded_data),], design = tumor_normal~ age_range + gender  + stages )

# Run DESeq2
dds <- DESeq(dds)

# Obtain differential expression numbers
resultsNames(dds)

# Save results for each stage
Stage_I   <-results(dds,name="stages_Stage.I_vs_")
Stage_II  <-results(dds,name="stages_Stage.II_vs_")
Stage_III <-results(dds,name="stages_Stage.III_vs_")

