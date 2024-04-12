#######################################################################################################################################
# Caclulate Gene expression pvalue and log2fc among samples of each group                                                             #
#######################################################################################################################################
colData_file                            <-"/home/felipe/Documentos/LungPortal/samples/colData.tsv"                                    #
unstranded_file                         <-"/home/felipe/Documentos/LungPortal/samples/unstranded_data_id.tsv"                         #
log2fc_expression_pos_file              <-"/home/felipe/Documentos/LungPortal/samples/log2fc_expression_pos.tsv"                      #
#######################################################################################################################################
colData_data                            <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)                        #
unstranded_data                         <-read.table(file = unstranded_file, sep = '\t', header = TRUE,fill=TRUE)                     #
log2fc_expression_pos_data              <-read.table(file = unstranded_file, sep = '\t', header = TRUE,fill=TRUE)                     #
#######################################################################################################################################
# Filter up genes                                                                                                                     #
unstranded_data<-unstranded_data[rownames(log2fc_expression_pos),colData_data$patient_id]                                             #
                                                                                                                                      #
# Filter up patientt                                                                                                                  #
colData_data<-colData_data[colnames(unstranded_data),]                                                                                #
                                                                                                                                      #
# Define all stage pairs                                                                                                              #
# Paired stages                                                                                                                       #                                                                                                                    #
stages_pairs=data.frame(stage_i=c("stageI","stageI","stageII"),stage_ii=c("stageII","stageIII","stageIII"))                           #
#######################################################################################################################################
# Genes of each stage stored in colData                                                                                               #
sample_stage_I  <-colData[colData$stages=="Stage I","patient_id"]                                                                     #
sample_stage_II <-colData[colData$stages=="Stage II","patient_id"]                                                                    #
sample_stage_III<-colData[colData$stages=="Stage III","patient_id"]                                                                   #
#######################################################################################################################################
# Lists for stage names,samples, and genes.                                                                                           #
vector_stages   <- c("stageI","stageII","stageIII")                                                                                   #
list_samples    <- list(stageI=sample_stage_I,stageII=sample_stage_II,stageIII=sample_stage_III)                                      #
#######################################################################################################################################

# A table for each gene, with the following columns:
# # 1-Gene  2-StageI_StageII_log2foldchange  3-StageI_StageII_pvalue 4-StageI_StageIII_log2foldchange  5-StageI_StageIII_pvalue  6-StageII_StageIII_log2foldchange  7-StageII_StageIII_pvalue
# log2foldchange = log2(expression value in condition A) - log2(expression value in condition B)

# For each genes, complete the t.test pvalue and log2foldchange for any pair of stages
for (gene in rownames(unstranded_data))
{

  # For each stage pair                                                                                                                                                                                                                                      #
  for (stage_pair in rownames(stages_pairs))                                                                                                                                                                                                                 #
  {                                                                                                                                                                                                                                                          #
      # Store pairs                                                                                                                                                                                                                                          #
      stage_i <- stages_pairs[stage_pair,"stage_i"]                                                                                                                                                                                                          #
      stage_ii<- stages_pairs[stage_pair,"stage_ii"]      

    

      unstranded_data[,]
  }
    
  
}
