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
unstranded_data               <-read.table(file = unstranded_data_file, sep = '\t', header = TRUE,fill=TRUE)         #
merged_data_patient_info_data      <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)#
colData_data                       <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)                 #
rownames(colData)                  <-colData$patient_id                                                                   #
###########################################################################################################################
dds_stages_tissue_type <- DESeqDataSetFromMatrix(countData = unstranded_data, colData=colData_data[colnames(unstranded_data),], design = ~  age_at_index +  gender +tissue_type  )
###########################################################################################################################

# Paired samples
paired_sample_df<-data.frame(normal=c(normal_samples_id),tumor=c(tumor_solid_sample_id),case=case)

# For each case, find the pairs
for (case in unique(merged_data_patient_info_data$case))
{
    # Take the samples of this case (normal and tumor)
    case_samples<-merged_data_patient_info_data[merged_data_patient_info_data$case==case,]

    # Take the tumor samples
    tumor_sampĺes <-case_samples[case_samples$tissue_type=="Tumor",]
    normal_sampĺes<-case_samples[case_samples$tissue_type=="Normal",]

    # if vector contains at least one tumor and one normal
    if(length(unique(normal_sampĺes$sample_id))>0 && length(unique(tumor_sampĺes$sample_id))>0)
    {
            # For each tumor sample
            for (tumor_solid_sample_id in tumor_sampĺes$sample_id)
            {
                # for each normal sample, compile a paired samples
                for (normal_samples_id in normal_sampĺes$sample_id)
                {
                    # Contatenate                     
                    paired_sample_df<-rbind(data.frame(normal=c(normal_samples_id),tumor=c(tumor_solid_sample_id),case=case),paired_sample_df)
                }
            }                
    }
}
