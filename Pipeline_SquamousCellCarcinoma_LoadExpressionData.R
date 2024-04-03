#####################################################################################################################
# This script will take the TSV file (metadata), unstranded.rna_seq.augmented_star_gene_counts (rna-seq count data), 
# the name of genes and samples already processed for the primary_diagnosis=Squamous cell carcinoma, NOS
#####################################################################################################################
library(readr)
library("xlsx")
library(ggplot2)
library("DESeq2")
library(gridExtra)
#####################################################################################################################
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

# Abril 03 00:31
# Here I created the field patient_id
# I called as such, the world patient is used
# but what I meant is, some cases have more than one sample.
# it can be what I call patient_id is one one entry for the sample (with the file)
# Check what it is this of samples_ids that cannot be used
# How : check the field /home/felipe/Documentos/LungPortal/samples/merged_data_patient_info.tsv
# I suspect the inconsistence is not incosistence
# Create field patient_id
merged_data_patient_info$patient_id<-paste("patient_",1:length(merged_data_patient_info$sample_id),sep="")

# Filter collumns that are used for age_at_index, gender, stages, Sample.ID
merged_data_patient_info<-merged_data_patient_info[,c("patient_id","case_id","sample_id","age_at_index","gender","Sample.Type","stages","race","tissue_type")]

# Rename collumns
colnames(merged_data_patient_info)<-c("patient_id","case_id","sample_id","age_at_index","gender","tumor_normal","stages","race","tissue_type")
#####################################################################################################################
# Set colData
colData<-unique(merged_data_patient_info[,c("sample_id","age_at_index","gender","tumor_normal","stages","race","tissue_type")])

# Patient and sample
patient_sample<-merged_data_patient_info[which(merged_data_patient_info$sample_id %in% colData$sample_id),c("patient_id","sample_id")]

# Take firtst occurance
patient_sample_first <- patient_sample[match(unique(patient_sample$sample_id), patient_sample$sample_id),]

patient_sample_first
colData$patient_id<-patient_sample_first[which(patient_sample_first$sample_id %in% colData$sample_id),"patient_id"]

# Set colnames
rownames(colData)<-colData$sample_id
#####################################################################################################################
# To DO:
# Filter merged_data_patient_info stages with N>30
merged_data_patient_info<-merged_data_patient_info[merged_data_patient_info$stage %in%  names(which(table(merged_data_patient_info$stage)>30)),]

# Filter merged_data_patient_info stages with N>30
unstranded_data<-unstranded_data[,colnames(unstranded_data) %in% merged_data_patient_info$sample_id]

# Filter colData
colData<-colData[unique(merged_data_patient_info$sample_id),]

# Filter DESeq2 steps
#####################################################################################################################
# Set colnames
rownames(colData)<-colData$sample_id

# Set rownames
colnames(unstranded_data) <- colData[colnames(unstranded_data),"patient_id"]

# Set colnames
rownames(colData)<-colData$patient_id
#####################################################################################################################
# https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2019.00930/full#h3
# The gene expression data were obtained as RNA-seq files in their version 2 (Illumina Hi-Seq) available for tissues affected by cancer or not (paired tissues), 
# from TCGA (https://cancergenome.nih.gov/) accessed in February 2016. Version 2 gives gene expression values for 20,532 genes referred to as GeneSymbol, 
# calculated by RNA-seq through expectation maximization (RSEM) (Li and Dewey, 2011) and normalized according to the upper quartile methods. 
# The 9,190 genes for which the equivalence between GeneSymbols and UniProtKB could be obtained went through further analysis. 
# This equivalence list is available in Supplementary Table 1. 
# To do : use only genes in the list
# /home/felipe/Documentos/LungPortal/Table_1.
gene_ids_file       <- "/home/felipe/Documentos/LungPortal/samples/gene_ids.txt"                            
gene_name_file      <- "/home/felipe/Documentos/LungPortal/samples/gene_name.txt"

# Load gene data (name and id)
gene_ids_data<-read.table(file = gene_ids_file, sep = '\t', header = FALSE,fill=TRUE)      
gene_name_data<-read.table(file = gene_name_file, sep = '\t', header = FALSE,fill=TRUE)      

# Put both in a data
df_gene_id_symbol<-data.frame(gene_id=gene_ids_data,gene_symbol=gene_name_data)

# Rename collumns
colnames(df_gene_id_symbol)<-c("gene_id","gene_symbol")

library("readxl")
Table1_data<-read.table(file = "/home/felipe/Documentos/LungPortal/Table_1.tsv", sep = '\t', header = TRUE,fill=TRUE)      

# Genes in conforte et al data
selected_gene_id<-df_gene_id_symbol[df_gene_id_symbol$gene_symbol %in% Table1_data$GeneSymbol,"gene_id"]
#####################################################################################################################
write.table(unstranded_data, file = "/home/felipe/Documentos/LungPortal/samples/unstranded.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
#####################################################################################################################
# Filter RNA-seq data to contain only data from Conforte et.al
unstranded_data<-unstranded_data[selected_gene_id,]
#####################################################################################################################
merged_data_patient_info$patient_id<-paste("patient_",1:length(merged_data_patient_info$sample_id),sep="")

# Writing mtcars data
write.table(colData, file = "/home/felipe/Documentos/LungPortal/samples/colData.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(merged_data_patient_info[,c("patient_id","case_id","sample_id","age_at_index","gender","tumor_normal","stages","tissue_type")], file ="/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv" , sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(unstranded_data, file = "/home/felipe/Documentos/LungPortal/samples/unstranded_data_id.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
#####################################################################################################################
