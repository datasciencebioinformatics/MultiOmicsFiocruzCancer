##########################################################################################################################################################################################################
library(readr)
library("xlsx")
library(ggplot2)
# A table containing the Transcriptome Profiling of each sample
# Read the sample tsv file and for each "Transcriptome Profiling" read the file
 # Reading the contents of TSV file using read_tsv() method
gdc_sample_sheet_file<-"/home/felipe/Documentos/LungPortal/gdc_sample_sheet.2024-03-08.tsv"

# Read data
gdc_sample_sheet_data<-read.table(file = gdc_sample_sheet_file, sep = '\t', header = TRUE,fill=TRUE)  

# Add collumn sample_id
gdc_sample_sheet_data$sample_submitter_id<-gdc_sample_sheet_data$Sample.ID
#####################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output/"                                                           #
#####################################################################################################################
#To do:
# - Merge data info with patient info
# - Associationg between Normal/Tumor and Covariables.
# - Check criterias for selectig data e.g. "Transcriptome Profiling" only. e.g. ".FPKM.txt.gz"
# - Outputs : a) a table dataInfo + patientInfo.
# -           b) A table with expression of selected patients.
# - To plan : 1) association covariable~primary_diagnosis
# -           2) association covariable~Normal/Tumor
# - To study : Literature of association co-variables ~ cancer
#####################################################################################################################
# Set path to files
clinical_file="/home/felipe/Documentos/LungPortal/clinical.txt" 
sample_file="/home/felipe/Documentos/LungPortal/sample.txt"    
exposure_file="/home/felipe/Documentos/LungPortal/exposure.txt"                                                  #

# Load data
clinical_data<-read.table(file = clinical_file, sep = '\t', header = TRUE,fill=TRUE)    
sample_data<-read.table(file = sample_file, sep = '\t', header = TRUE,fill=TRUE)                                    #
exposure_data<-read.table(file = exposure_file, sep = '\t', header = TRUE,fill=TRUE)                                #

# Merge data
merged_sample_clinical_data<-merge(sample_data,clinical_data,by="case_id")

# Merge all
merged_sample_clinical_data<-merge(merged_sample_clinical_data,exposure_data,by="case_id")

# Merge tables
merged_data_patient_info<-merge(merged_sample_clinical_data,gdc_sample_sheet_data,by="sample_submitter_id")
#####################################################################################################################
# Set file name variable 
merged_data_patient_info$<-merged_data_patient_info$File.Name


#####################################################################################################################
# A script to load cancer data base in R
unstranded_file       <- "/home/felipe/Documentos/LungPortal/samples/unstranded.rna_seq.augmented_star_gene_counts.tsv"
colnames_file         <- "/home/felipe/Documentos/LungPortal/samples/header.txt"
rownames_file         <- "/home/felipe/Documentos/LungPortal/samples/gene_ids.txt"

# Load data
unstranded_data<-read.table(file = unstranded_file, sep = '\t', header = FALSE,fill=TRUE)    
colnames_data<-read.table(file = colnames_file, sep = '\t', header = FALSE,fill=TRUE)                                    #
rownames_data<-read.table(file = rownames_file, sep = '\t', header = FALSE,fill=TRUE)      

# Set colnames and rownames
rownames(unstranded_data)<-rownames_data[,1]
colnames(unstranded_data)<-colnames_data[,1]

# To check counts and sample_ids
sum(merged_data_patient_info$sample_id %in% colnames(unstranded_data))
