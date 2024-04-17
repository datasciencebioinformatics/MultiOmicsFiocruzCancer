library("reshape2")
library("ggVennDiagram")
##########################################################################################################################################################
# A scrtpt to calculate                                                                                                                                  #
# 1) boxplots for DE genes from each stage (e.g. Stage I vs. Steges II and III) using smples of each Stage (e.g. Samples from Stage I),                  #
# and also from the genes in the intersection of the stages                                                                                              #
# 2) order of magnitudes among any two stages norm_counts(stage_ii)/norm_counts(stage_i)                                                                 #
# 3) correlation among order of magniture and number of connections from interactome                                                                     #
# obs: one positive genes selected from both DE normal/tumor and DE stages are used                                                                      #
#                                                                                                                                                        #
##########################################################################################################################################################
# Path to files of selected_genes                                                                                                      # 
sele_genes_uniq_pos_stages_I_file       <-"/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_pos_stage_I.tsv"             #
sele_genes_uniq_pos_stages_II_file      <-"/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_pos_stage_II.tsv"            #
sele_genes_uniq_pos_stages_III_file     <-"/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_pos_stage_III.tsv"           #
unstranded_file                         <-"/home/felipe/Documentos/LungPortal/samples/unstranded_data_id.tsv"                          #
colData_file                            <-"/home/felipe/Documentos/LungPortal/samples/colData.tsv"                                     #
df_log2foldchange_file                  <-"/home/felipe/Documentos/LungPortal/samples/df_log2foldchange.tsv"                           #
########################################################################################################################################
# Load data                                                                                                                           #
sele_genes_uniq_pos_stages_I_data       <-read.table(file = sele_genes_uniq_pos_stages_I_file, sep = '\t', header = TRUE,fill=TRUE)   #
sele_genes_uniq_pos_stages_II_data      <-read.table(file = sele_genes_uniq_pos_stages_II_file, sep = '\t', header = TRUE,fill=TRUE)  #
sele_genes_uniq_pos_stages_III_data     <-read.table(file = sele_genes_uniq_pos_stages_III_file, sep = '\t', header = TRUE,fill=TRUE) #
unstranded_data                         <-read.table(file = unstranded_file, sep = '\t', header = TRUE,fill=TRUE)                     #
colData_data                            <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)                        #
df_log2foldchange_data                  <-read.table(file = df_log2foldchange_file, sep = '\t', header = TRUE,fill=TRUE)              #
                                                                                                                                      #
# Set rownames                                                                                                                        #
rownames(sele_genes_uniq_pos_stages_I_data)<-sele_genes_uniq_pos_stages_I_data$Gene                                                   #
rownames(sele_genes_uniq_pos_stages_II_data)<-sele_genes_uniq_pos_stages_II_data$Gene                                                 #
rownames(sele_genes_uniq_pos_stages_III_data)<-sele_genes_uniq_pos_stages_III_data$Gene                                               #
#######################################################################################################################################
# Select stages 
stages_I_II_III<-ggVennDiagram(list(Stage_I    =sele_genes_uniq_pos_stages_I_data$Gene,Stage_II  =sele_genes_uniq_pos_stages_II_data$Gene,Stage_III  =sele_genes_uniq_pos_stages_III_data$Gene), label_alpha = 0) + scale_fill_viridis() + theme_bw() + ggtitle("Stages I, II and III")
#######################################################################################################################################
# Filter to keep only positive tumor/normal samples from unstranded_data.                                                             #
unstranded_data<-unstranded_data[avg_expression_pos$Gene,]                                                                            #
                                                                                                                                      #
# omit NA values from vector                                                                                                          #
unstranded_data <- na.omit(unstranded_data)                                                                                           #
                                                                                                                                      #
# Filter to keep only positive tumor/normal samples from colData_data.                                                                #
colData_data<-colData_data[colnames(unstranded_data),]                                                                                #
                                                                                                                                      #                                                                                                                                      #
# Create dds_stages with selected for the purpose of getting normalized counts                                                        #
dds_stages <- DESeqDataSetFromMatrix(countData = unstranded_data, colData=colData, design = ~  age_at_index +  gender +tissue_type  ) #
                                                                                                                                      #
# Estimate size factor                                                                                                                #
dds_stages <- estimateSizeFactors(dds_stages)                                                                                         #
                                                                                                                                      #
# Obtain normalized coutns                                                                                                            #
norm_counts<-counts(dds_stages, normalized = TRUE)                                                                                    #
#######################################################################################################################################
# Genes of each stage stored in colData                                                                                               #
sample_stage_I  <-colData[colData$stages=="Stage I","patient_id"]                                                                     #
sample_stage_II <-colData[colData$stages=="Stage II","patient_id"]                                                                    #
sample_stage_III<-colData[colData$stages=="Stage III","patient_id"]                                                                   #
#######################################################################################################################################
merged_data_patient_info_file      <- "/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv"                               #
merged_data_patient_info_data      <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)            #
########################################################################################################################################################
# Lists for stage names,samples, and genes.                                                                                                            #
vector_stages   <- c("stageI","stageII","stageIII")                                                                                                    #
list_samples    <- list(stageI=sample_stage_I,stageII=sample_stage_II,stageIII=sample_stage_III)                                                       #
list_genes      <- list(stageI=sele_genes_uniq_pos_stages_I_data$Gene,stageII=sele_genes_uniq_pos_stages_II_data$Gene,stageIII=sele_genes_uniq_pos_stages_III_data$Gene) #
############################################################################################################################################################################################################################################################
# Paired stages                                                                                                                                                                                                                                            #
stages_pairs=data.frame(stage_i=c("stageI","stageI","stageII"),stage_ii=c("stageII","stageIII","stageIII"))                                                                                                                                                #
                                                                                                                                                                                                                                                           #
# Data frame results                                                                                                                                                                                                                                       #
df_order_of_magnitude<-data.frame(Genes=c(),stage_pair=c(),order_of_magnitude=c())                                                                                                                                                                         #
                                                                                                                                                                                                                                                           #
# For each stage pair                                                                                                                                                                                                                                      #
for (stage_pair in rownames(stages_pairs))                                                                                                                                                                                                                 #
{                                                                                                                                                                                                                                                          #
    # Store pairs                                                                                                                                                                                                                                          #
    stage_i <- stages_pairs[stage_pair,"stage_i"]                                                                                                                                                                                                          #
    stage_ii<- stages_pairs[stage_pair,"stage_ii"]                                                                                                                                                                                                         #
                                                                                                                                                                                                                                                           #
    # Store pairs                                                                                                                                                                                                                                          #
    genes_stage_i  <- list_genes[[stage_i]]                                                                                                                                                                                                                #
    genes_stage_ii <- list_genes[[stage_ii]]                                                                                                                                                                                                               #
                                                                                                                                                                                                                                                           # 
    # Take the samples and genes of stages                                                                                                                                                                                                                 #
    samples_stage_i  <-list_samples[[stage_i]]                                                                                                                                                                                                             #
    samples_stage_ii <-list_samples[[stage_ii]]       
                                                                                                                                                           #
    # Paste Star pairs                                                                                                                                                                                                                                     ############################################################
    df_order_of_magnitude<-rbind(df_order_of_magnitude,data.frame(Genes=names(rowMeans(norm_counts[,samples_stage_ii])/rowMeans(norm_counts[,samples_stage_i])),stage_pair=paste(stage_ii,stage_i,sep="_over_"),order_of_magnitude=rowMeans(norm_counts[,samples_stage_ii])/rowMeans(norm_counts[,samples_stage_i]))) #
}                                                                                                                                                                                                                                                          ############################################################
############################################################################################################################################################################################################################################################
library("ggpubr")
# Create log2 order of magnitude                                                                                                                                                                                                                           #
df_order_of_magnitude$log2_order_of_magnitude<-log(df_order_of_magnitude$order_of_magnitude,2)                                                                                                                                                             #
                                                                                                                                                                                                                                                           #
# Pairwise comparisons: Specify the comparisons you want                                                                                                                                                                                                   #
my_comparisons <- list( c("stageII_over_stageI", "stageIII_over_stageI"), c("stageIII_over_stageII", "stageIII_over_stageI"), c("stageIII_over_stageII", "stageII_over_stageI") )                                                                          #                                                                                                                                                                                                                                                         #
                                                                                                                                                                                                                                                           #
############################################################################################################################################################################################################################################################
# df_order_of_magnitude                                                                                              #
# Selcted only gennes that can converted in the table of PPI                                                         #
df_order_of_magnitude<-df_order_of_magnitude[df_order_of_magnitude$Gene %in% merge_interactome_gene_symbol$gene_id,]                   #
                                                                                                                     #                                                                                                                   #
#for each rown, add the number of connections                                                                        #
df_order_of_magnitude$PPI<-0                                                                                         #

# Subset gene_symbol and PPI
df_interactome<-unique(data.frame(Gene=merge_interactome_gene_symbol$gene_id,PPI=merge_interactome_gene_symbol$PPI))

# Take just the first occurance
df_interactome<-df_interactome[!duplicated(df_interactome$Gene), ]

# Set rownames()
rownames(df_interactome)<-df_interactome$Gene

# Assert PPI
df_order_of_magnitude$PPI<-df_interactome[df_order_of_magnitude$Genes,"PPI"]                                         

df_magnitude_stageI<-cbind(df_order_of_magnitude,data.frame(Stage="Stage I",log2foldchange=df_log2foldchange_data[as.vector(df_order_of_magnitude$Genes),"StageI_StagesII_III_log2foldchange2"]))
df_magnitude_stageII<-cbind(df_order_of_magnitude,data.frame(Stage="Stage II",log2foldchange=df_log2foldchange_data[as.vector(df_order_of_magnitude$Genes),"StageII_StagesI_III_log2foldchange2"]))
df_magnitude_stageIII<-cbind(df_order_of_magnitude,data.frame(Stage="Stage III",log2foldchange=df_log2foldchange_data[as.vector(df_order_of_magnitude$Genes),"StageIII_StagesI_II_log2foldchange"]))
                                                                                               
# Assert Average Gene Expression
df_order_of_magnitude_melt<-rbind(df_magnitude_stageI,df_magnitude_stageII,df_magnitude_stageIII)
                                                                                                                   #
# PPI                                                                                                                #
df_order_of_magnitude_melt<-na.omit(df_order_of_magnitude_melt)                                                                #

cor(df_order_of_magnitude_melt$Epxr_Mean, df_order_of_magnitude_melt$PPI, method = c("pearson", "kendall", "spearman"))
cor(df_magnitude_stageI$Epxr_Mean, df_magnitude_stageI$PPI, method = c("pearson", "kendall", "spearman"))
cor(df_magnitude_stageII$Epxr_Mean, df_magnitude_stageII$PPI, method = c("pearson", "kendall", "spearman"))
cor(df_magnitude_stageIII$Epxr_Mean, df_magnitude_stageIII$PPI, method = c("pearson", "kendall", "spearman"))


                                                                                                                     #
# FindClusters_resolution                                                                                            #
png(filename=paste(output_dir,"log2foldchange_vs_Order_Magnitude.png",sep=""), width = 20, height = 14, res=600, units = "cm")     ########################################################################################################
  ggplot(df_order_of_magnitude_melt, aes(y =log(order_of_magnitude), x = log(log2foldchange))) +  geom_point() + theme_bw()  + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) + ggtitle("log2foldchange vs. Order of magnitude")   + facet_wrap(~stage_pair, ncol = 3)         
dev.off()   

# FindClusters_resolution                                                                                            #
png(filename=paste(output_dir,"PPI_vs_log2foldchange.png",sep=""), width = 20, height = 14, res=600, units = "cm")     ########################################################################################################
  ggplot(df_order_of_magnitude_melt, aes(x = log(PPI), y = log2foldchange)) +  geom_point() + theme_bw()  + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) + ggtitle("Connectivity vs. log2foldchange")   + facet_wrap(~Stage, ncol = 3)          +stat_cor(method = "pearson") 
dev.off()   
