# A R script to create metadata from gdc files.
# Inputs:
# gdc_sample_sheet.2024-03-08.tsv
# /home/felipe/Documentos/LungPortal/clinical.txt
# /home/felipe/Documentos/LungPortal/sample.txt
# /home/felipe/Documentos/LungPortal/exposure.txt
# Output : merged_data_patient_info.tsv
###########################################################################################################################
merged_data_patient_info_file       <- "/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv"                  #
colData_file                        <- "/home/felipe/Documentos/LungPortal/samples/colData.tsv"                           #
###########################################################################################################################
merged_data_patient_info_data      <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)#
colData_data                       <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)                 #
rownames(colData)                  <-colData$patient_id                                                                   #
#######################################################################################################################################
colData<-na.omit(colData)                                                                                                             #
merged_data_patient_info_data<-merged_data_patient_info_data[which(merged_data_patient_info_data$patient_id %in% colData$patient_id),]#
#######################################################################################################################################
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
#######################################################################################################################################
# folchange=Expr(Stage i)/Expr(Stage ii and II)
#folchange=rowMeans(unstranded_rpkm[,paired_sample_df$tumor])-rowMeans(unstranded_rpkm[,paired_sample_df$normal])
folchange=rowMeans(unstranded_rpkm[,colData[colData$tissue_type=="Tumor","patient_id"]])/rowMeans(unstranded_rpkm[,paired_sample_df$normal])

# log2change
log2change=log(folchange,2)	

# log2change data
log2change_tumor_control=na.omit(data.frame(gene=names(log2change),log2change=log2change))
#######################################################################################################################################
