library(readr)
library("xlsx")
library(ggplot2)
##########################################################################################################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output/"                                                             #
#####################################################################################################################
# Reading the contents of TSV file using read_tsv() method
merged_data_patient_info_file<- "/home/felipe/Documentos/LungPortal/samples/merged_data_patient_info.tsv"

# Load metadata table
merged_data_patient_info     <-read.table(file = merged_data_patient_info, sep = '\t', header = TRUE,fill=TRUE)   
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
# There are 238 cases, 570 samples, 478 tumor, 92 normal

# Filter collumns that are used for age_at_index, gender, stages, Sample.ID
merged_data_patient_info<-merged_data_patient_info[,c("sample_id","age_at_index","gender","stages")]

## and to read this file back into R one needs
read.table("/home/felipe/Documentos/LungPortal/samples/unstranded.rna_seq.gene_counts.tsv", header = TRUE, sep = ",", row.names = 1)

# Organize how to send to Carles
write_tsv(unstranded_data, "/home/felipe/Documentos/LungPortal/samples/unstranded.rna_seq.gene_counts.tsv")
write_tsv(merged_data_patient_info, "/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv")

# Creat DEseq element from unstranded_data and merged_data_patient_info
dds <- DESeqDataSetFromMatrix(countData = unstranded_data, colData=merged_data_patient_info, design = ~sample_id + age_at_index + gender + stages)

# Run DESeq2
dds <- DESeq(dds)

# Obtain differential expression numbers
res <- results(dds, alpha=0.05, lfcThreshold=log2(2))
