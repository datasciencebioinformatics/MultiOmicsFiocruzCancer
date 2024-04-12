########################################################################################################################################
# Caclulate Gene expression pvalue and log2fc among samples of each group                                                              #
########################################################################################################################################
colData_file                            <-"/home/felipe/Documentos/LungPortal/samples/colData.tsv"                                     #
unstranded_file                         <-"/home/felipe/Documentos/LungPortal/samples/unstranded_data_id.tsv"                          #
log2fc_expression_pos_file              <-"/home/felipe/Documentos/LungPortal/samples/log2fc_expression_pos.tsv"                       #
########################################################################################################################################
colData_data                            <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)                         #
unstranded_data                         <-read.table(file = unstranded_file, sep = '\t', header = TRUE,fill=TRUE)                      #
log2fc_expression_pos_data              <-read.table(file = unstranded_file, sep = '\t', header = TRUE,fill=TRUE)                      #
########################################################################################################################################
# Filter up genes                                                                                                                      #
unstranded_data<-unstranded_data[rownames(log2fc_expression_pos),colData_data$patient_id]                                              #
                                                                                                                                       #
# Filter up patientt                                                                                                                   #
colData_data<-colData_data[colnames(unstranded_data),]                                                                                 #
########################################################################################################################################

# Define all stage pairs 

# Take the samples of each stage ()

# For each genes, complete the t.test pvalue and log2foldchange for any pair of stages
