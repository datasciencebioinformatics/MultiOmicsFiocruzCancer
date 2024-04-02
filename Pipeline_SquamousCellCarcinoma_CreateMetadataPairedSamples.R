# A R script to create metadata from gdc files.
# Inputs:
# gdc_sample_sheet.2024-03-08.tsv
# /home/felipe/Documentos/LungPortal/clinical.txt
# /home/felipe/Documentos/LungPortal/sample.txt
# /home/felipe/Documentos/LungPortal/exposure.txt
# Output : merged_data_patient_info.tsv
##########################################################################################################################################################################################################
library(readr)
library("xlsx")
library(ggplot2)
##########################################################################################################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output/"                                                                   #
###########################################################################################################################
unstranded_data_file                <- "/home/felipe/Documentos/LungPortal/samples/unstranded_data_id.tsv"                #
merged_data_patient_info_file       <- "/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv"                  #
colData_file                        <- "/home/felipe/Documentos/LungPortal/samples/colData.tsv"                           #
###########################################################################################################################
unstranded_data_data               <-read.table(file = unstranded_data_file, sep = '\t', header = TRUE,fill=TRUE)         #
merged_data_patient_info_data      <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)#
colData_data                       <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)                 #
rownames(colData)                  <-colData$patient_id                                                                   #
###########################################################################################################################
# For each case, find the pairs
for (case in unique(merged_data_patient_info_data$case))
{
    # Take the samples of this case (normal and tumor)
    case_samples<-merged_data_patient_info_data[merged_data_patient_info_data$case==case,]

    # Take the tumor samples
    tumor_sampĺes <-case_samples[case_samples$tissue_type=="Tumor",]
    normal_sampĺes<-case_samples[case_samples$tissue_type=="Normal",]

    # For each tumor sample
    for (tumor_sample in unique(tumor_sampĺes$sample_id))
    {
        
    }        
}
