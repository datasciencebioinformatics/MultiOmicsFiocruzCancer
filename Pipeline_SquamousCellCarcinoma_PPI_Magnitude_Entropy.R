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
########################################################################################################################################
# Load data                                                                                                                           #
sele_genes_uniq_pos_stages_I_data       <-read.table(file = sele_genes_uniq_pos_stages_I_file, sep = '\t', header = TRUE,fill=TRUE)   #
sele_genes_uniq_pos_stages_II_data      <-read.table(file = sele_genes_uniq_pos_stages_II_file, sep = '\t', header = TRUE,fill=TRUE)  #
sele_genes_uniq_pos_stages_III_data     <-read.table(file = sele_genes_uniq_pos_stages_III_file, sep = '\t', header = TRUE,fill=TRUE) #
unstranded_data                         <-read.table(file = unstranded_file, sep = '\t', header = TRUE,fill=TRUE)                     #
colData_data                            <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)                        #
                                                                                                                                      #
# Set rownames                                                                                                                        #
rownames(sele_genes_uniq_pos_stages_I_data)<-sele_genes_uniq_pos_stages_I_data$Gene                                                   #
rownames(sele_genes_uniq_pos_stages_II_data)<-sele_genes_uniq_pos_stages_II_data$Gene                                                 #
rownames(sele_genes_uniq_pos_stages_III_data)<-sele_genes_uniq_pos_stages_III_data$Gene                                               #
#######################################################################################################################################
intersection_genes_pos_Stages_I_II_III_file <- "/home/felipe/Documentos/LungPortal/output/intersection_genes_pos_Stages_I_II_III.tsv" #
intersection_genes_pos_Stages_I_II_file     <- "/home/felipe/Documentos/LungPortal/output/intersection_genes_pos_Stages_I_II.tsv"     #
intersection_genes_pos_Stages_I_III_file    <- "/home/felipe/Documentos/LungPortal/output/intersection_genes_pos_Stages_I_III.tsv"    #
intersection_genes_pos_Stages_II_III_file   <- "/home/felipe/Documentos/LungPortal/output/intersection_genes_pos_Stages_II_III.tsv"   #
######################################################################################################################################################
# Load data                                                                                                                                          #
intersection_genes_pos_Stages_I_II_III       <-read.table(file = intersection_genes_pos_Stages_I_II_III_file, sep = '\t', header = TRUE,fill=TRUE)   #
intersection_genes_pos_Stages_I_II           <-read.table(file = intersection_genes_pos_Stages_I_II_file, sep = '\t', header = TRUE,fill=TRUE)       #
intersection_genes_pos_Stages_II_III         <-read.table(file = intersection_genes_pos_Stages_II_III_file, sep = '\t', header = TRUE,fill=TRUE)     #
intersection_genes_pos_Stages_I_III          <-read.table(file = intersection_genes_pos_Stages_I_III_file, sep = '\t', header = TRUE,fill=TRUE)      #
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
                                                                                                                                      #
                                                                                                                                      #
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
########################################################################################################################################################################################################################
# A table with DE genes from each stage (e.g. Stage I vs. Steges II and III) using smples of each Stage (e.g. Samples from Stage I)                                                                                    #                                                                                                                     #
# DE Genes Stage I                                                                                                                                                                                                     #
DE_genes_from_stage_I_vs_Samples_from_stage_I   <- data.frame(DE_genes_stage="DE genes from stage I",samples="Samples from stage I", Expr=melt(norm_counts[sele_genes_uniq_pos_stages_I_data$Gene,sample_stage_I]))          #
DE_genes_from_stage_I_vs_Samples_from_stage_II  <- data.frame(DE_genes_stage="DE genes from stage I",samples="Samples from stage II", Expr=melt(norm_counts[sele_genes_uniq_pos_stages_I_data$Gene,sample_stage_II]))        #
DE_genes_from_stage_I_vs_Samples_from_stage_III <- data.frame(DE_genes_stage="DE genes from stage I",samples="Samples from stage III", Expr=melt(norm_counts[sele_genes_uniq_pos_stages_I_data$Gene,sample_stage_III]))      #
                                                                                                                                                                                                                       #
# DE Genes Stage II                                                                                                                                                                                                    #
DE_genes_from_stage_II_vs_Samples_from_stage_I   <- data.frame(DE_genes_stage="DE genes from stage II",samples="Samples from stage I", Expr=melt(norm_counts[sele_genes_uniq_pos_stages_II_data$Gene,sample_stage_I]))       #
DE_genes_from_stage_II_vs_Samples_from_stage_II  <- data.frame(DE_genes_stage="DE genes from stage II",samples="Samples from stage II", Expr=melt(norm_counts[sele_genes_uniq_pos_stages_II_data$Gene,sample_stage_II]))     #
DE_genes_from_stage_II_vs_Samples_from_stage_III <- data.frame(DE_genes_stage="DE genes from stage II",samples="Samples from stage III", Expr=melt(norm_counts[sele_genes_uniq_pos_stages_II_data$Gene,sample_stage_III]))   #
                                                                                                                                                                                                                       #
# DE Genes Stage III                                                                                                                                                                                                   #
DE_genes_from_stage_III_vs_Samples_from_stage_I   <- data.frame(DE_genes_stage="DE genes from stage III",samples="Samples from stage I", Expr=melt(norm_counts[sele_genes_uniq_pos_stages_III_data$Gene,sample_stage_I]))    #
DE_genes_from_stage_III_vs_Samples_from_stage_II  <- data.frame(DE_genes_stage="DE genes from stage III",samples="Samples from stage II", Expr=melt(norm_counts[sele_genes_uniq_pos_stages_III_data$Gene,sample_stage_II]))  #
DE_genes_from_stage_III_vs_Samples_from_stage_III <- data.frame(DE_genes_stage="DE genes from stage III",samples="Samples from stage III", Expr=melt(norm_counts[sele_genes_uniq_pos_stages_III_data$Gene,sample_stage_III]))#
                                                                                                                                                                                                                       #
# Combine rbind tables in single table                                                                                                                                                                                 #
melt_combined_table<-rbind(DE_genes_from_stage_I_vs_Samples_from_stage_I,DE_genes_from_stage_I_vs_Samples_from_stage_II,DE_genes_from_stage_I_vs_Samples_from_stage_III,DE_genes_from_stage_II_vs_Samples_from_stage_I,
DE_genes_from_stage_II_vs_Samples_from_stage_II, DE_genes_from_stage_II_vs_Samples_from_stage_III, DE_genes_from_stage_III_vs_Samples_from_stage_I,DE_genes_from_stage_III_vs_Samples_from_stage_II,                  DE_genes_from_stage_III_vs_Samples_from_stage_III )                                                                                                                                                                    #
                                                                                                                                                                                                                       #
# Create log Expr.value                                                                                                                                                                                                #
melt_combined_table$log2Expr<-log(melt_combined_table$Expr.value,2)                                                                                                                                                    #
                                                                                                                                                                                                                       #
# Rename collumns                                                                                                                                                                                                      #
colnames(melt_combined_table)<-c("DE_genes_stage", "samples", "Genes", "Patients", "log2Expr")                                                                                                                         #
                                                                                                                                                                                                                       #
# plot                                                                                                                                                                                                                 #
p1 <- ggplot(melt_combined_table, aes(x=DE_genes_stage, y=log2Expr, fill=samples)) +                                                                                                                                   #
    geom_boxplot(varwidth = TRUE) + facet_grid(~samples) +                                                                                                                                                             #
    theme(legend.position="none") + geom_boxplot() + scale_fill_brewer(palette="Set1") + theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +                                                 #
    ggtitle("Expression of genes per stage") + ylim(0, 10)                                                                                                                                                             #
                                                                                                                                                                                                                       # 
# FindClusters_resolution                                                                                                                                                                                              #
png(filename=paste(output_dir,"Expression_of_genes_per_stage.png",sep=""), width = 18, height = 18, res=600, units = "cm")                                                                                             #
  p1 
dev.off()                                                                                                                                                                                                              #
############################################################################################################################################################################################################################################################
# Samples from stage I                                                                                                                                                                                                                                     #
DE_genes_from_Intersection_I_II_III_vs_Samples_from_stage_I   <-  data.frame(DE_genes_stage="Intersection stages I, II, III",samples="Samples from stage I", Expr=melt(norm_counts[intersection_genes_pos_Stages_I_II_III$Gene,sample_stage_I]))           #
DE_genes_from_Intersection_I_II_vs_Samples_from_stage_I       <-  data.frame(DE_genes_stage="Intersection stages I, II",samples="Samples from stage I", Expr=melt(norm_counts[intersection_genes_pos_Stages_I_II_III$Gene,sample_stage_I]))                #
DE_genes_from_Intersection_II_III_vs_Samples_from_stage_I     <-  data.frame(DE_genes_stage="Intersection stages II, III",samples="Samples from stage I", Expr=melt(norm_counts[intersection_genes_pos_Stages_I_II_III$Gene,sample_stage_I]))              #
                                                                                                                                                                                                                                                           #
# Samples from stage II                                                                                                                                                                                                                                    #
DE_genes_from_Intersection_I_II_III_vs_Samples_from_stage_II   <-  data.frame(DE_genes_stage="Intersection stages I, II, III",samples="Samples from stage II", Expr=melt(norm_counts[intersection_genes_pos_Stages_I_II_III$Gene,sample_stage_II]))        #
DE_genes_from_Intersection_I_II_vs_Samples_from_stage_II       <-  data.frame(DE_genes_stage="Intersection stages I, II",samples="Samples from stage II", Expr=melt(norm_counts[intersection_genes_pos_Stages_I_II_III$Gene,sample_stage_II]))             #
DE_genes_from_Intersection_II_III_vs_Samples_from_stage_II     <-  data.frame(DE_genes_stage="Intersection stages II, III",samples="Samples from stage II", Expr=melt(norm_counts[intersection_genes_pos_Stages_I_II_III$Gene,sample_stage_II]))           #
                                                                                                                                                                                                                                                           #
# Samples from stage III                                                                                                                                                                                                                                   #
DE_genes_from_Intersection_I_II_III_vs_Samples_from_stage_III  <-  data.frame(DE_genes_stage="Intersection stages I, II, III",samples="Samples from stage III", Expr=melt(norm_counts[intersection_genes_pos_Stages_I_II_III$Gene,sample_stage_III]))      #
DE_genes_from_Intersection_I_II_vs_Samples_from_stage_III      <-  data.frame(DE_genes_stage="Intersection stages I, II",samples="Samples from stage III", Expr=melt(norm_counts[intersection_genes_pos_Stages_I_II_III$Gene,sample_stage_III]))           #
DE_genes_from_Intersection_II_II_vs_Samples_from_stage_III     <-  data.frame(DE_genes_stage="Intersection stages II, III",samples="Samples from stage III", Expr=melt(norm_counts[intersection_genes_pos_Stages_I_II_III$Gene,sample_stage_III]))         #
                                                                                                                                                                                                                                                           #
# Combine rbind tables in single table                                                                                                                                                                                                                     #
melt_combined_table_intersection<-rbind(DE_genes_from_Intersection_I_II_III_vs_Samples_from_stage_I,DE_genes_from_Intersection_I_II_vs_Samples_from_stage_I,DE_genes_from_Intersection_II_III_vs_Samples_from_stage_I,                                     #
DE_genes_from_Intersection_I_II_III_vs_Samples_from_stage_II,DE_genes_from_Intersection_I_II_vs_Samples_from_stage_II,DE_genes_from_Intersection_II_III_vs_Samples_from_stage_II,DE_genes_from_Intersection_I_II_III_vs_Samples_from_stage_III,            #
DE_genes_from_Intersection_I_II_vs_Samples_from_stage_III,DE_genes_from_Intersection_II_II_vs_Samples_from_stage_III)                                                                                                                                      #
                                                                                                                                                                                                                                                           #
# Add combined table intersection                                                                                                                                                                                                                          #
melt_combined_table_intersection$log2Expr<-log(melt_combined_table_intersection$Expr.value)                                                                                                                                                                #
                                                                                                                                                                                                                                                           #
# plot                                                                                                                                                                                                                                                     #
p2 <- ggplot(melt_combined_table_intersection, aes(x=DE_genes_stage, y=log2Expr, fill=samples)) +                                                                                                                                                          #
    geom_boxplot(varwidth = TRUE) + facet_grid(~samples) +                                                                                                                                                                                                 #
    theme(legend.position="none") + geom_boxplot() + scale_fill_brewer(palette="Set1") + theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) + ggtitle("Expression of genes per stage")                                            #
############################################################################################################################################################################################################################################################                                                                                                     #
# FindClusters_resolution                                                                                                                                                                                                                                  #
png(filename=paste(output_dir,"Expression_of_genes_per_intersaction.png",sep=""), width = 18, height = 18, res=600, units = "cm")                                                                                                                          #
  p2                                                                                                                                                                                                                                                       #
dev.off()                                                                                                                                                                                                                                                  #
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
    samples_stage_ii <-list_samples[[stage_ii]]                                                                                                                                                                                                            #
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

# Set rownames()
rownames(df_interactome)<-df_interactome$Gene

# Take just the first occurance
df_interactome<-df_interactome[!duplicated(df_interactome$Gene), ]

# Assert rownames
rownames(df_interactome)<-df_interactome$Gene

# Assert PPI
df_order_of_magnitude$PPI<-df_interactome[df_order_of_magnitude$Genes,"PPI"]                                         

# Assert Average Gene Expression
df_order_of_magnitude$Epxr_Mean<-rowMeans(norm_counts[df_order_of_magnitude$Genes,])
                                                                                                                     #
# PPI                                                                                                                #
df_order_of_magnitude<-na.omit(df_order_of_magnitude)                                                                #
                                                                                                                     #
# FindClusters_resolution                                                                                            #
png(filename=paste(output_dir,"PPI_vs_Order_Magnitude.png",sep=""), width = 20, height = 14, res=600, units = "cm")     ########################################################################################################
  ggplot(df_order_of_magnitude, aes(x = log(PPI), y = log(order_of_magnitude))) +  geom_point() + theme_bw()  + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) + ggtitle("PPI vs. Order of magnitude")   + facet_wrap(~stage_pair, ncol = 3)          +stat_cor(method = "pearson")  # Add correlation coefficient
dev.off()   

# FindClusters_resolution                                                                                            #
png(filename=paste(output_dir,"PPI_vs_Expression.png",sep=""), width = 20, height = 14, res=600, units = "cm")     ########################################################################################################
  ggplot(df_order_of_magnitude, aes(x = log(PPI), y = log(Epxr_Mean))) +  geom_point() + theme_bw()  + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) + ggtitle("PPI vs. Epxr_Mean")   + facet_wrap(~stage_pair, ncol = 3)          +stat_cor(method = "pearson") 
dev.off()   

#############################################################################################################################################################################################################################                                                                                                      
# Set rownames                                                                                                                                                                                                              #
combined_genes_intersection<-unique(c(rownames(sele_genes_uniq_pos_stages_I_data),rownames(sele_genes_uniq_pos_stages_II_data),rownames(sele_genes_uniq_pos_stages_III_data)))                                              #
                                                                                                                                                                                                                            #
# Combine all intersections                                                                                                                                                                                                 #
combined_genes_intersection<-unique(c(intersection_genes_pos_Stages_I_II_III$Gene,intersection_genes_pos_Stages_I_II$Gene,intersection_genes_pos_Stages_II_III$Gene,intersection_genes_pos_Stages_I_III$Gene))             #
                                                                                                                                                                                                                            #
# Comb_genes_from_Intersections_vs_Samples_from_stage                                                                                                                                                                       #
Comb_genes_from_Intersections_vs_Samples_from_stage_I     <-  data.frame(DE_genes_stage="Intersection stages I, II, III",stages="Stage I", Expr=melt(norm_counts[combined_genes_intersection,sample_stage_I]))              #
Comb_genes_from_Intersections_vs_Samples_from_stage_II    <-  data.frame(DE_genes_stage="Intersection stages I, II, III",stages="Stage II", Expr=melt(norm_counts[combined_genes_intersection,sample_stage_II]))            #
Comb_genes_from_Intersections_vs_Samples_from_stage_III   <-  data.frame(DE_genes_stage="Intersection stages I, II, III",stages="Stage III", Expr=melt(norm_counts[combined_genes_intersection,sample_stage_III]))          #
                                                                                                                                                                                                                            #
# melt_comb_genes                                                                                                                                                                                                           #
melt_comb_genes_from_Intersections<-rbind(Comb_genes_from_Intersections_vs_Samples_from_stage_I,Comb_genes_from_Intersections_vs_Samples_from_stage_II,Comb_genes_from_Intersections_vs_Samples_from_stage_III)             #
                                                                                                                                                                                                                            #
                                                                                                                                                                                                                            #
# Create log Expr.value                                                                                                                                                                                                     #
melt_comb_genes_from_Intersections$log2Expr<-log(melt_comb_genes_from_Intersections$Expr.value,2)                                                                                                                           #                               #
                                                                                                                                                                                                                            #
# Rename collumns                                                                                                                                                                                                           #
colnames(melt_comb_genes_from_Intersections)<-c("DE_genes_stage", "stages", "Genes", "Patients", "log2Expr")                                                                                                                #
                                                                                                                                                                                                                            #
# plot                                                                                                                                                                                                                      #
p4 <- ggplot(melt_comb_genes_from_Intersections, aes(x=stages, y=log2Expr, fill=stages)) + geom_boxplot(varwidth = TRUE) +  theme(legend.position="none") + geom_boxplot() + scale_fill_brewer(palette="Set1") + theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) + ggtitle("Expression of DE genes from all genes combined")    + ylim(0, 5)    
                                                                                                                                                                                                                            #
# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"Combined_dE_genes_per_stage.png",sep=""), width = 18, height = 18, res=600, units = "cm")                                                                                                    #
  p4                                                                                                                                                                                                                        #
dev.off()                                                                                                                                                                                                                   #
#########################################################################################################################################################################################
# To do : to check samples in the intesctons of the stage I would have to cck for patient ids with samples samples fom more than one stage                                              #                    
# To do : to select top 10 genes and plot paired overview, boxplot and individually.                                                                                                    # 
#########################################################################################################################################################################################
# merged_data_patient_info_sub
merged_data_patient_info_sub<-merged_data_patient_info_data[,-(1)]

merged_data_patient_info_sub<-

cases_id_stageI<-unique(merged_data_patient_info_sub[merged_data_patient_info_sub$stages=="Stage I","case_id"])
cases_id_stageII<-unique(merged_data_patient_info_sub[merged_data_patient_info_sub$stages=="Stage II","case_id"])
cases_id_stageIII<-unique(merged_data_patient_info_sub[merged_data_patient_info_sub$stages=="Stage III","case_id"])

# Combine gnes from each estage
df_StageI_Samples<-rbind(data.frame(Expr=melt(norm_counts[intersection_genes_pos_Stages_I_II_III$Gene,sample_stage_I]),Stage="Stage I", Samples="Stage I"),
data.frame(Expr=melt(norm_counts[intersection_genes_pos_Stages_I_II_III$Gene,sample_stage_II]),Stage="Stage I", Samples="Stage II"),
data.frame(Expr=melt(norm_counts[intersection_genes_pos_Stages_I_II_III$Gene,sample_stage_III]),Stage="Stage I", Samples="Stage III"))

# Rename collumns                                                                                                                                                                                                           #
colnames(df_StageI_Samples)<-c("Gene", "Patient", "Expr", "Stage","Samples")  

# Organize only top genes
df_StageI_Samples$log2Expr<-log(df_StageI_Samples$Expr)

# Basic barplot
p5<-ggplot(data=df_StageI_Samples, aes(x=Samples, y=log2Expr)) +  geom_bar(stat="identity") + facet_wrap(~Gene, ncol = 3) + scale_fill_brewer(palette="Set1") + theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) + ggtitle("DE genes for selected genes present in samples from stage I, II and III")   

# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"melt_comb_genes_from_Intersection_I_II_II.png",sep=""), width = 18, height = 18, res=600, units = "cm")                                                                                                    #
  p5
dev.off() 
#########################################################################################################################################################################################
# Pairwise comparisons: Specify the comparisons you want                                                                                                                                                                                                   #
my_comparisons <- list( c("Stage I", "Stage II"), c("Stage II", "Stage III"),c("Stage I", "Stage III"))                                                                                 #

# ggboxplot
p <- ggboxplot(df_StageI_Samples, x = "Samples", y = "log2Expr", color = "Samples", palette = "jco",  add = "jitter")  + stat_compare_means(comparisons = my_comparisons)   + scale_fill_brewer(palette="Set1") + theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) + ggtitle("DE genes for selected genes present in samples from stage I, II and III")   + facet_wrap(~Gene, ncol = 3, scales = "free") 

# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"melt_comb_genes_intersaction_pvalus.png",sep=""), width = 24, height = 24, res=600, units = "cm")                                                                                                    #
  p
dev.off() 
