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

# Data.frame df_log2foldchange
df_log2foldchange<-data.frame(Gene=rownames(unstranded_data),StageI_StageII_log2foldchange=0,StageI_StageII_pvalue=0,StageI_StageIII_log2foldchange=0,StageI_StageIII_pvalue=0,StageII_StageIII_log2foldchange=0,StageII_StageIII_pvalue=0)

# Set rownames
rownames(df_log2foldchange)<-df_log2foldchange$Gene

# For each genes, complete the t.test pvalue and log2foldchange for any pair of stages
for (gene in rownames(unstranded_data))
{

  # index_log2foldchange and index_pvalue
  index_log2foldchange   <- 2
  index_pvalue           <- 3

  # For each stage pair                                                                                                                                                                                                                                      #
  for (stage_pair in rownames(stages_pairs))                                                                                                                                                                                                                 #
  {                                                                                                                                                                                                                                                          #
      # Store pairs                                                                                                                                                                                                                                          #
      stage_i  <- list_samples[stages_pairs[stage_pair,"stage_i"]]
      stage_ii <- list_samples[stages_pairs[stage_pair,"stage_ii"]]

      # Take gene expresion in each group
      gene_expression_stage_i <-unstranded_data[gene,as.vector(unlist(stage_i))]
      gene_expression_stage_ii<-unstranded_data[gene,as.vector(unlist(stage_ii))]

      # Calculate log2foldchange
      #calulate the average values in each group
      mean_stage_i  = rowMeans(gene_expression_stage_i)
      mean_stage_ii = rowMeans(gene_expression_stage_ii)

      # logFC and pvalue
      logFC=log2(mean_stage_i/mean_stage_ii)
      pvalue=t.test(gene_expression_stage_i, gene_expression_stage_ii, alternative = "two.sided", var.equal = FALSE)$p.value

      df_log2foldchange[gene,index_log2foldchange]<-logFC
      df_log2foldchange[gene,index_pvalue]        <-pvalue
      df_log2foldchange[gene,"Gene"]              <-gene  

      # index_log2foldchange and index_pvalue
      index_log2foldchange   <- index_log2foldchange+2
      index_pvalue           <- index_pvalue+2    
  }  
}
