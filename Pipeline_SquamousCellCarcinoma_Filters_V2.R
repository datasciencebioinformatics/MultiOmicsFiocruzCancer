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

# Save tumor samples
tumor_samples<-colData[which(colData$tissue_type=="Tumor"),"patient_id"]

# Filter by log2folchage
unstranded_rpkm<-unstranded_rpkm[rownames(log2change_tumor_control[which(log2change_tumor_control$log2change>=log2foldchange_trheshold),]),]

# Filter by RPKM
unstranded_data_filter<-unstranded_rpkm[rowMeans(unstranded_rpkm[,tumor_samples])>10,]
###########################################################################################################################
print(paste("Number of up-regulated tumor-genes :",dim(unstranded_data_filter)[1]))
cat(print(paste("\nNumber of up-regulated tumor-genes :",dim(unstranded_data_filter)[1])),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)

