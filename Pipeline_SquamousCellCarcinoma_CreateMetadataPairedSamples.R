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
unstranded_data_file                <- "/home/felipe/Documentos/LungPortal/samples/unstranded.tsv"                        #
merged_data_patient_info_file       <- "/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv"                  #
colData_file                        <- "/home/felipe/Documentos/LungPortal/samples/colData.tsv"                           #
###########################################################################################################################
unstranded_data                    <-read.table(file = unstranded_data_file, sep = '\t', header = TRUE,fill=TRUE)         #
merged_data_patient_info_data      <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)#
colData_data                       <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)                 #
rownames(colData)                  <-colData$patient_id                 
############################################################################################################################
colData<-na.omit(colData)
unstranded_data<-unstranded_data[,colData$patient_id]
merged_data_patient_info_data<-merged_data_patient_info_data[which(merged_data_patient_info_data$patient_id %in% colData$patient_id),]
###########################################################################################################################
dds_stages_tissue_type <- DESeqDataSetFromMatrix(countData = unstranded_data, colData=colData_data[colnames(unstranded_data),], design = ~  age_at_index +  gender +tissue_type  )
###########################################################################################################################
# Estimate size factor
dds_stages_tissue_type <- estimateSizeFactors(dds_stages_tissue_type)

# Obtain normalized coutns
norm_counts<-counts(dds_stages_tissue_type, normalized = TRUE)
###########################################################################################################################
# Paired samples
paired_sample_df<-data.frame(normal=c(),tumor=c(),case=c())

# For each case, find the pairs
for (case in unique(merged_data_patient_info_data$case))
{
    # All samples for case id = "case"
    case_samples<-merged_data_patient_info_data[merged_data_patient_info_data$case==case,]

    # Take the tumor samples
    tumor_sampĺes <-case_samples[case_samples$tissue_type=="Tumor",]
    normal_sampĺes<-case_samples[case_samples$tissue_type=="Normal",]

    # if vector contains at least one tumor and one normal
    if(length(unique(normal_sampĺes$sample_id))>0 && length(unique(tumor_sampĺes$sample_id))>0)
    {
            # For each tumor sample
            for (tumor_solid_sample_id in tumor_sampĺes$patient_id)
            {
                # for each normal sample, compile a paired samples
                for (normal_samples_id in normal_sampĺes$patient_id)
                {
                    # Contatenate                     
                    paired_sample_df<-rbind(data.frame(normal=c(normal_samples_id),tumor=c(tumor_solid_sample_id),case=case),paired_sample_df)
                }
            }                
    }
}
###########################################################################################################################
# Data.frame to store results for control samples
df_diff_expression<-data.frame(Gene=rownames(unstranded_data))

# set rownames
rownames(df_diff_expression)<-df_diff_expression$Gene

# For each pair, we subtracted gene expression values of control samples
for (case in rownames(paired_sample_df))
{    

    case_sample  <-paired_sample_df[case,"case"]
    normal_sample<-paired_sample_df[case,"normal"]
    tumor_sample <-paired_sample_df[case,"tumor"]    
   
    # Normal and tumor samples for this "case" id
    normal_sample<-merged_data_patient_info_data[which(merged_data_patient_info_data$patient_id==normal_sample),]
    tumor_sample <-merged_data_patient_info_data[which(merged_data_patient_info_data$patient_id==tumor_sample),]

    # Store results for the case
    case_results_normal<-data.frame(case=norm_counts[,c(normal_sample$patient_id)])
    case_results_tumor <-data.frame(case=norm_counts[,c(tumor_sample$patient_id)])
    diff_expression<-case_results_tumor-case_results_normal

    # Rename collumns
    colnames(diff_expression)<-case_sample

    # Concatenate 
    df_diff_expression<-cbind(df_diff_expression,diff_expression)    
}
# Remove gene expression
df_diff_expression<-df_diff_expression[,-1]

# Take average expression
avg_expression<-data.frame(rowMeans(df_diff_expression))

# Take average expression of positive sample
avg_expression_pos<-data.frame(Gene=rownames(avg_expression)[which(avg_expression>=0.5)],Expression=avg_expression[which(avg_expression>=0.5),])

# Take average expression of positive sample
rownames(avg_expression_pos)<-avg_expression_pos$Gene

# Write table
write.table(avg_expression_pos, file = "/home/felipe/Documentos/LungPortal/samples/avg_expression_pos.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)


