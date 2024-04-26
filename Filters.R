###########################################################################################################################
unstranded_data_file                <- "/home/felipe/Documentos/LungPortal/samples/unstranded_rpkm.tsv"                    #
merged_data_patient_info_file       <- "/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv"                  #
colData_file                        <- "/home/felipe/Documentos/LungPortal/samples/colData.tsv"                           #
###########################################################################################################################
unstranded_data                    <-read.table(file = unstranded_data_file, sep = '\t', header = TRUE,fill=TRUE)         #
merged_data_patient_info_data      <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)#
colData_data                       <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)                 #
rownames(colData)                  <-colData$patient_id                                                                   #
###########################################################################################################################
tumor_samples<-colData[colData$tissue_type=="Tumor","patient_id"]
normal_samples<-colData[colData$tissue_type=="Normal","patient_id"]

# Filter by rowmeans(tumor)-rowmeans(normal)
unstranded_data_filter<-unstranded_data[rowMeans(unstranded_data[,tumor_samples])-rowMeans(unstranded_data[,normal_samples])>0,]

# Filter by rowmeans(tumor)-rowmeans(normal)
unstranded_data_filter<-unstranded_data[rowMeans(unstranded_data[,tumor_samples])>10,]
###########################################################################################################################
# Save TSV file with genes from Stage1                                                                                                                             #
write.table(data.frame(unstranded_rpkm), file =  "/home/felipe/Documentos/LungPortal/samples/unstranded_rpkm.tsv", append = FALSE,row.names = TRUE, col.names = TRUE, quote = TRUE, sep = "\t")
###########################################################################################################################


