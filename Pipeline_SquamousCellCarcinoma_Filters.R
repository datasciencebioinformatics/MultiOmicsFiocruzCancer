###########################################################################################################################
merged_data_patient_info_file       <- "/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv"                  #
colData_file                        <- "/home/felipe/Documentos/LungPortal/samples/colData.tsv"                           #
###########################################################################################################################
unstranded_data                    <- unstranded_rpkm
merged_data_patient_info_data      <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)#
colData_data                       <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)                 #
rownames(colData)                  <-colData$patient_id                                                                   #
###########################################################################################################################
log2foldchange_trheshold<-log2foldchange_trheshold
# log2fc_threshold
log2fc_threshold<-log2fc_threshold

tumor_samples<-colData_data[colData_data$tissue_type=="Tumor","patient_id"]
normal_samples<-colData_data[colData_data$tissue_type=="Normal","patient_id"]

# Filter by rowmeans(tumor)-rowmeans(normal)
log2fc<-log(rowMeans(unstranded_rpkm[,tumor_samples]),2)/log(rowMeans(unstranded_rpkm[,normal_samples],2))

# Filter by log2folchage
unstranded_rpkm<-unstranded_rpkm[log2fc>log2foldchange_trheshold,]

# Filter by RPKM
unstranded_data_filter<-unstranded_rpkm[rowMeans(unstranded_rpkm[,tumor_samples])>10,]
###########################################################################################################################
print(paste("Number of up-regulated tumor-genes :",dim(unstranded_data_filter)[1]))

