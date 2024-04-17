######################################################################################################################################
# Path to files of selected_genes                                                                                                    #
unstranded_file                         <-"/home/felipe/Documentos/LungPortal/samples/unstranded_data_id.tsv"                        #
colData_file                            <-"/home/felipe/Documentos/LungPortal/samples/colData.tsv"                                   #
######################################################################################################################################
# Load data                                                                                                                          #
unstranded_data                         <-read.table(file = unstranded_file, sep = '\t', header = TRUE,fill=TRUE)                    #
colData_data                            <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)                       #
######################################################################################################################################
# Genes of each stage stored in colData                                                                                              #
sample_stage_I  <-colData[colData$stages=="Stage I","patient_id"]                                                                    #
sample_stage_II <-colData[colData$stages=="Stage II","patient_id"]                                                                   #
sample_stage_III<-colData[colData$stages=="Stage III","patient_id"]                                                                  #
######################################################################################################################################
# Genes of each stage stored in colData                                                                                              #
sample_stage_I  <-colData[colData$stages=="Stage I","patient_id"]                                                                    #
sample_stage_II  <-colData[colData$stages=="Stage II","patient_id"]                                                                  #
sample_stage_III  <-colData[colData$stages=="Stage III","patient_id"]                                                                #
######################################################################################################################################
colData_stage_I  <-colData[sample_stage_I,]                                                                                          #
colData_stage_II <-colData[sample_stage_II,]                                                                                         #
colData_stage_III<-colData[sample_stage_III,]                                                                                        #
                                                                                                                                     #
rawdata_stage_I  <-unstranded_data[,sample_stage_I]                                                                                  #
rawdata_stage_II <-unstranded_data[,sample_stage_II]                                                                                 #
rawdata_stage_III<-unstranded_data[,sample_stage_III]                                                                                #
###################################################################################################################################################################
# Write table                                                                                                                                                     #
write.table(colData_stage_I, file = "/home/felipe/Documentos/LungPortal/output/stages/colData_stage_I.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)       #
write.table(colData_stage_II, file = "/home/felipe/Documentos/LungPortal/output/stages/colData_stage_II.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)     #
write.table(colData_stage_III, file = "/home/felipe/Documentos/LungPortal/output/stages/colData_stage_III.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)   #
                                                                                                                                                                  #
# Write table                                                                                                                                                     #
write.table(rawdata_stage_I, file = "/home/felipe/Documentos/LungPortal/output/stages/rawdata_stage_I.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)       #
write.table(rawdata_stage_II, file = "/home/felipe/Documentos/LungPortal/output/stages/rawdata_stage_II.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)     #
write.table(rawdata_stage_III, file = "/home/felipe/Documentos/LungPortal/output/stages/rawdata_stage_III.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)   #
###################################################################################################################################################################
# Re-Load data                                                                                                                       #
unstranded_data                         <-read.table(file = unstranded_file, sep = '\t', header = TRUE,fill=TRUE)                    #
colData_data                            <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)                       #
######################################################################################################################################
# Separate in normal and tumor raw counts                                               #
normal_raw_counts  <-unstranded_data[,paired_sample_df$normal]                          #
normal_tumor_counts<-unstranded_data[,paired_sample_df$tumor]                           #
                                                                                        #
# Subselect colData                                                                     #
colData_sub<-colData[c(colnames(normal_raw_counts),colnames(normal_tumor_counts)),]     #
                                                                                        #
# Rename collumns                                                                       #
colnames(normal_raw_counts)  <-paired_sample_df$case                                    #
colnames(normal_tumor_counts)<-paired_sample_df$case                                    #
                                                                                        #
# Set cases                                                                             # 
colData_sub[paired_sample_df$normal,"case"]<-paired_sample_df$case                      #
colData_sub[paired_sample_df$tumor,"case"] <-paired_sample_df$case                      # 
                                                                                        #
# Separate normal and tumor                                                             #
colData_normal<-colData_sub[paired_sample_df$normal,]                                   #
colData_tumor<-colData_sub[paired_sample_df$tumor,]                                     #
                                                                                        #
# Genes of each stage stored in colData                                                 #   
sample_stage_I_normal  <-colData_normal[colData_normal$stages=="Stage I","patient_id"]  #
sample_stage_II_normal <-colData_normal[colData_normal$stages=="Stage II","patient_id"] #
sample_stage_III_normal<-colData_normal[colData_normal$stages=="Stage III","patient_id"]#
                                                                                        #
# Genes of each stage stored in colData                                                 #
sample_stage_I_tumor  <-colData_tumor[colData_tumor$stages=="Stage I","patient_id"]     #
sample_stage_II_tumor <-colData_tumor[colData_tumor$stages=="Stage II","patient_id"]    #
sample_stage_III_tumor<-colData_tumor[colData_tumor$stages=="Stage III","patient_id"]   #    
#########################################################################################
colData_normal_stage_I<-colData_normal[sample_stage_I_normal,]                         #
colData_normal_stage_II<-colData_normal[sample_stage_II_normal,]                       #
colData_normal_stage_III<-colData_normal[sample_stage_III_normal,]                     #
                                                                                       #
colData_tumor_stage_I<-colData_tumor[sample_stage_I_tumor,]                            #
colData_tumor_stage_II<-colData_tumor[sample_stage_II_tumor,]                          #
colData_tumor_stage_III<-colData_tumor[sample_stage_III_tumor,]                        #
################################################################################################
rawdata_normal_stage_I<-unstranded_data[,sample_stage_I_normal]                                #
rawdata_normal_stage_II<-unstranded_data[,sample_stage_II_normal]                              #
rawdata_normal_stage_III<-unstranded_data[,sample_stage_III_normal]                            #
                                                                                               #
colnames(rawdata_normal_stage_I)   <- colData_normal[colnames(rawdata_normal_stage_I),"case"]  # 
colnames(rawdata_normal_stage_II)  <- colData_normal[colnames(rawdata_normal_stage_II),"case"] #
colnames(rawdata_normal_stage_III) <- colData_normal[colnames(rawdata_normal_stage_III),"case"]#
                                                                                               #
rawdata_tumor_stage_I<-unstranded_data[,sample_stage_I_tumor]                                  #
rawdata_tumor_stage_II<-unstranded_data[,sample_stage_II_tumor]                                #
rawdata_tumor_stage_III<-unstranded_data[,sample_stage_III_tumor]                              #
                                                                                               #
colnames(rawdata_tumor_stage_I)   <- colData_tumor[colnames(rawdata_tumor_stage_I),"case"]     #
colnames(rawdata_tumor_stage_II)  <- colData_tumor[colnames(rawdata_tumor_stage_II),"case"]    #
colnames(rawdata_tumor_stage_III) <- colData_tumor[colnames(rawdata_tumor_stage_III),"case"]   #
####################################################################################################################################################################################
# Write table                                                                                                                                                                      #
write.table(colData_normal_stage_I, file = "/home/felipe/Documentos/LungPortal/output/stages/colData_normal_stage_I.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)          #
write.table(colData_normal_stage_II, file = "/home/felipe/Documentos/LungPortal/output/stages/colData_normal_stage_II.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)        #
write.table(colData_normal_stage_III, file = "/home/felipe/Documentos/LungPortal/output/stages/colData_normal_stage_III.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)      #
                                                                                                                                                                                   #
# Write table                                                                                                                                                                      #
write.table(colData_tumor_stage_I, file = "/home/felipe/Documentos/LungPortal/output/stages/colData_tumor_stage_I.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)            #
write.table(colData_tumor_stage_II, file = "/home/felipe/Documentos/LungPortal/output/stages/colData_tumor_stage_II.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)          #
write.table(colData_tumor_stage_III, file = "/home/felipe/Documentos/LungPortal/output/stages/colData_tumor_stage_III.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)        #
                                                                                                                                                                                   #
write.table(rawdata_normal_stage_I, file = "/home/felipe/Documentos/LungPortal/output/stages/rawdata_normal_stage_I.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)          #
write.table(rawdata_normal_stage_II, file = "/home/felipe/Documentos/LungPortal/output/stages/rawdata_normal_stage_II.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)        #
write.table(rawdata_normal_stage_III, file = "/home/felipe/Documentos/LungPortal/output/stages/rawdata_normal_stage_III.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)      #
                                                                                                                                                                                   #
write.table(rawdata_tumor_stage_I, file = "/home/felipe/Documentos/LungPortal/output/stages/rawdata_tumor_stage_I.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)            #
write.table(rawdata_tumor_stage_II, file = "/home/felipe/Documentos/LungPortal/output/stages/rawdata_tumor_stage_II.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)          #
write.table(rawdata_tumor_stage_III, file = "/home/felipe/Documentos/LungPortal/output/stages/rawdata_tumor_stage_III.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)        #
####################################################################################################################################################################################
