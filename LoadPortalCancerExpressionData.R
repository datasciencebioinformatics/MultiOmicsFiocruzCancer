library(readr)
library("xlsx")
library(ggplot2)
library("DESeq2")
library(gridExtra)
#####################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output/"                                                             #
#####################################################################################################################
# Reading the contents of TSV file using read_tsv() method
merged_data_patient_info_file<- "/home/felipe/Documentos/LungPortal/samples/merged_data_patient_info.tsv"

# Load metadata table
merged_data_patient_info     <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)   
#####################################################################################################################
# A script to load cancer data base in R
unstranded_file       <- "/home/felipe/Documentos/LungPortal/samples/unstranded.rna_seq.augmented_star_gene_counts.tsv"
colnames_file         <- "/home/felipe/Documentos/LungPortal/samples/header.txt"
rownames_file         <- "/home/felipe/Documentos/LungPortal/samples/gene_ids.txt"


# Load data
unstranded_data<-read.table(file = unstranded_file, sep = '\t', header = FALSE,fill=TRUE)    
colnames_data<-read.table(file = colnames_file, sep = '\t', header = FALSE,fill=TRUE)                                    
rownames_data<-read.table(file = rownames_file, sep = '\t', header = FALSE,fill=TRUE)      

# Set colnames and rownames
rownames(unstranded_data)<-rownames_data[,1]
colnames(unstranded_data)<-colnames_data[,1]

# To check counts and sample_ids
sum(merged_data_patient_info$sample_id %in% colnames(unstranded_data))

# A field to store 
merged_data_patient_info$stages<-merged_data_patient_info$ajcc_pathologic_stage

# Group stages I,II,III and IV
merged_data_patient_info$stages<-gsub("Stage IA", "Stage I", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IB", "Stage I", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIA", "Stage II", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIB", "Stage II", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIIA", "Stage III", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIIB", "Stage III", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IVA", "Stage IV", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IVB", "Stage IV", merged_data_patient_info$stages)

# Filter the table to contain only samples in the dataset
merged_data_patient_info<-merged_data_patient_info[merged_data_patient_info$sample_id %in% colnames(unstranded_data),]

# Save table with samples and metatada
# There are 238 cases, 285 samples, 570 pairs

# Filter collumns that are used for age_at_index, gender, stages, Sample.ID
merged_data_patient_info<-merged_data_patient_info[,c("case_id","sample_id","age_at_index","gender","Sample.Type","stages")]

# Rename collumns
colnames(merged_data_patient_info)<-c("case_id","sample_id","age_at_index","gender","tumor_normal","stages")
#####################################################################################################################
# Set colData
colData<-unique(merged_data_patient_info[,c("sample_id","age_at_index","gender","tumor_normal","stages")])

# Set colnames
rownames(colData)<-colData$sample_id
#####################################################################################################################
# To DO:
# Filter merged_data_patient_info stages with N>30
merged_data_patient_info<-merged_data_patient_info[merged_data_patient_info$stage %in%  names(which(table(merged_data_patient_info$stage)>30)),]

# Filter merged_data_patient_info stages with N>30
unstranded_data<-unstranded_data[,colnames(unstranded_data) %in% merged_data_patient_info$sample_id]

# Filter colData
colData<-colData[unique(merged_data_patient_info$sample_id),]

# Filter DESeq2 steps
#####################################################################################################################
# Add columns for patient ID
colData$patient_id<-paste("patient_",1:length(colData$sample_id),sep="")

# Set colnames
rownames(colData)<-colData$sample_id

# Set rownames
colnames(unstranded_data) <- colData[colnames(unstranded_data),"patient_id"]

# Set colnames
rownames(colData)<-colData$patient_id
#####################################################################################################################
# https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2019.00930/full#h3
# The gene expression data were obtained as RNA-seq files in their version 2 (Illumina Hi-Seq) available for tissues affected by cancer or not (paired tissues), 
# from TCGA (https://cancergenome.nih.gov/) accessed in February 2016. Version 2 gives gene expression values for 20,532 genes referred to as GeneSymbol, 
# calculated by RNA-seq through expectation maximization (RSEM) (Li and Dewey, 2011) and normalized according to the upper quartile methods. 
# The 9,190 genes for which the equivalence between GeneSymbols and UniProtKB could be obtained went through further analysis. 
# This equivalence list is available in Supplementary Table 1. 
# To do : use only genes in the list
# /home/felipe/Documentos/LungPortal/Table_1.tsv
gene_ids_file       <- "/home/felipe/Documentos/LungPortal/samples/gene_ids.txt"                            
gene_name_file      <- "/home/felipe/Documentos/LungPortal/samples/gene_name.txt"

# Load gene data (name and id)
gene_ids_data<-read.table(file = gene_ids_file, sep = '\t', header = FALSE,fill=TRUE)      
gene_name_data<-read.table(file = gene_name_file, sep = '\t', header = FALSE,fill=TRUE)      

# Put both in a data
df_gene_id_symbol<-data.frame(gene_id=gene_ids_data,gene_symbol=gene_name_data)

# Rename collumns
colnames(df_gene_id_symbol)<-c("gene_id","gene_symbol")

library("readxl")
Table1_data<-read.table(file = "/home/felipe/Documentos/LungPortal/Table_1.tsv", sep = '\t', header = TRUE,fill=TRUE)      

# Genes in conforte et al data
selected_gene_id<-df_gene_id_symbol[df_gene_id_symbol$gene_symbol %in% Table1_data$GeneSymbol,"gene_id"]

#####################################################################################################################
# Filter RNA-seq data to contain only data from Conforte et.al
unstranded_data<-unstranded_data[selected_gene_id,]

#####################################################################################################################
merged_data_patient_info$patient_id<-paste("patient_",1:length(merged_data_patient_info$sample_id),sep="")
write_tsv(unstranded_data, "/home/felipe/Documentos/LungPortal/samples/unstranded.rna_seq.gene_counts.tsv")
write_tsv(merged_data_patient_info[,c("patient_id","case_id","sample_id","age_at_index","gender","tumor_normal","stages")], "/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv")
#####################################################################################################################
# Set collumn to factor
#####################################################################################################################
# Set collumn to factor
colData$sample_id<-factor(colData$sample_id)
colData$tumor_normal<-factor(colData$tumor_normal)
colData$stages<-factor(colData$stages)
colData$gender<-factor(colData$gender)
colData$age_range<-factor(cut(colData$age_at_index, breaks = c(0, 25, 50,75,100 )))
colData$Tumor_Stage<-paste(colData$tumor_normal,colData$stages,sep="")

# Creat DEseq element from unstranded_data and merged_data_patient_info
# Here, tumor_normal is the outcome
# age_range + gender + stages are the covariables
# for one reason, DESeq is setting the last variable to group vs. all
# Two ways to check this information 1) check group vs. all, check to pca's

# Now, I do have the expression
dds <- DESeqDataSetFromMatrix(countData = unstranded_data, colData=colData[colnames(unstranded_data),], design = ~0 + stages:tumor_normal + stages + tumor_normal + age_range + gender )
dds <- DESeqDataSetFromMatrix(countData = unstranded_data, colData=colData[colnames(unstranded_data),], design = ~0 + stages + tumor_normal + age_range + gender )

# If you use a design of ~0 + covariable, then you will have a coefficient for each level of condition. 
# This is one of the cases where it helps to not have condition at the end of the design (for convenience, we often recommend to put condition at the end, but not in this case).
#####################################################################################################################
# Run DESeq2
dds <- DESeq(dds)

# Obtain differential expression numbers
resultsNames(dds)
#####################################################################################################################
# PCA analysis for dds
# DESIGN : design = ~0 + stages + tumor_normal + age_range + gender
# DATA   : expression (unstranded_data), metadata (colData)
#####################################################################################################################
### Transform counts for data visualization
vst_dds <- vst(dds, blind = TRUE, nsub = 1000, fitType = "parametric")

# Analysis of PCA 
## All Genes
### Plot PCA 
pca_tumor_normal<-plotPCA(vst_dds, intgroup="tumor_normal") + theme_bw() + ggtitle("tumor_normal")
pca_gender      <-plotPCA(vst_dds, intgroup="gender") + theme_bw()             + ggtitle("gender")
pca_age_range   <-plotPCA(vst_dds, intgroup="age_range")+ theme_bw()        + ggtitle("age_range")
pca_stages      <-plotPCA(vst_dds, intgroup="stages")+ theme_bw()              + ggtitle("stages")
pca_tumor_stages<-plotPCA(vst_dds, intgroup="Tumor_Stage")+ theme_bw()         + ggtitle("Stages + Primary diagnosis")

library(gridExtra)
library(cowplot)
pca_plots<-grid.arrange(pca_tumor_normal, pca_gender,pca_age_range,pca_stages,  nrow = 2)

# FindClusters_resolution
png(filename=paste(output_dir,"pca_plots.png",sep=""), width = 24, height = 36, res=600, units = "cm")
	plot_grid(pca_tumor_normal, pca_gender,pca_age_range,pca_stages, pca_tumor_stages,         ncol = 2, nrow = 3)
dev.off()
#####################################################################################################################
# Analysis of differential expression by stages, "one group vs. all others".
# Obtain differential expression numbers
Normal_Tumor<-data.frame(results(dds,name="tumor_normalSolid.Tissue.Normal"))

Stage_I    <-data.frame(results(dds,contrast=list(c("stagesStage.I"), c("stagesStage.II","stagesStage.III"))))
Stage_II   <-data.frame(results(dds,contrast=list(c("stagesStage.II"), c("stagesStage.I","stagesStage.III")))) 
Stage_III   <-data.frame(results(dds,contrast=list(c("stagesStage.III"), c("stagesStage.I","stagesStage.II")))) 

# Filter NA values
Stage_I_sub<-na.omit(Stage_I)
Stage_II_sub<-na.omit(Stage_II)
Stage_III_sub<-na.omit(Stage_III)
########################################################################################################################
log2fc_threshold  = 6
padj_treshold     = 0.01

Stage_I_sub   <-Stage_I_sub[which(abs(Stage_I_sub$log2FoldChange)>log2fc_threshold),]
Stage_II_sub  <-Stage_II_sub[which(abs(Stage_II_sub$log2FoldChange)>log2fc_threshold),]
Stage_III_sub <-Stage_III_sub[which(abs(Stage_III_sub$log2FoldChange)>log2fc_threshold),]

Stage_I_sub   <-Stage_I_sub[Stage_I_sub$padj<padj_treshold,]
Stage_II_sub  <-Stage_II_sub[Stage_II_sub$padj<padj_treshold,]
Stage_III_sub <-Stage_III_sub[Stage_III_sub$padj<padj_treshold,]

Stage_I_bck   <-Stage_I_sub
Stage_II_bck  <-Stage_II_sub
Stage_III_bck <-Stage_III_sub
########################################################################################################################
vst_Stage_I_sub<-varianceStabilizingTransformation(dds[rownames(Stage_I_sub),], blind = TRUE, fitType = "parametric")
vst_Stage_II_sub<-varianceStabilizingTransformation(dds[rownames(Stage_II_sub),], blind = TRUE, fitType = "parametric")
vst_Stage_III_sub<-varianceStabilizingTransformation(dds[rownames(Stage_III_sub),], blind = TRUE, fitType = "parametric")

pca_stageI<-plotPCA(vst_Stage_I_sub, intgroup="Tumor_Stage") + theme_bw() + ggtitle("Tumor_Stage : DE Genes from Stage I")
pca_stageII<-plotPCA(vst_Stage_II_sub, intgroup="Tumor_Stage") + theme_bw() + ggtitle("Tumor_Stage : DE Genes from Stage II")
pca_stageIII<-plotPCA(vst_Stage_III_sub, intgroup="Tumor_Stage") + theme_bw() + ggtitle("Tumor_Stage : DE Genes from Stage III")

# FindClusters_resolution
png(filename=paste(output_dir,"pca_stages_normal_cancer.png",sep=""), width = 36, height = 48, res=600, units = "cm")
	plot_grid(pca_stageI, pca_stageII,pca_stageIII, ncol = 2, nrow = 3)
dev.off()
#####################################################################################################################
# Save differential expression table
write_tsv(Stage_I_sub, "/home/felipe/Documentos/LungPortal/samples/stages_StageI.tsv")
write_tsv(Stage_II_sub, "/home/felipe/Documentos/LungPortal/samples/stages_StageII.tsv")
write_tsv(Stage_III_sub, "/home/felipe/Documentos/LungPortal/samples/stages_StageIII.tsv")
#####################################################################################################################
# Save differential expression table
write_tsv(unstranded_data[rownames(Stage_I_sub),], "/home/felipe/Documentos/LungPortal/samples/rnaseq_raw_counts_StageI.tsv")
write_tsv(unstranded_data[rownames(Stage_II_sub),], "/home/felipe/Documentos/LungPortal/samples/rnaseq_raw_counts_StageII.tsv")
write_tsv(unstranded_data[rownames(Stage_III_sub),], "/home/felipe/Documentos/LungPortal/samples/rnaseq_raw_counts_StageIII.tsv")
#####################################################################################################################
# Assert group for stage I samples in the form of Stage I vs. All Others
colData(vst_Stage_I_sub)$StageI_one_against_All<-factor(colData(vst_Stage_I_sub)$stages=="Stage I")
colData(vst_Stage_II_sub)$StageII_one_against_All<-factor(colData(vst_Stage_II_sub)$stages=="Stage II")
colData(vst_Stage_III_sub)$StageIII_one_against_All<-factor(colData(vst_Stage_III_sub)$stages=="Stage III")
#####################################################################################################################
Stage_I_sub   <-Stage_I_sub_bck
Stage_II_sub  <-Stage_II_sub_bck
Stage_III_sub <-Stage_III_sub_bck

colData(vst_Stage_I_sub)$StageI_one_against_All    <-""
colData(vst_Stage_II_sub)$StageII_one_against_All  <-""
colData(vst_Stage_III_sub)$StageIII_one_against_All<-""

colData(vst_Stage_I_sub)$StageI_one_against_All[which(colData(vst_Stage_I_sub)$stages=="Stage I")]<-"Stage I"
colData(vst_Stage_I_sub)$StageI_one_against_All[which(colData(vst_Stage_I_sub)$stages!="Stage I")]<-"All others"
colData(vst_Stage_II_sub)$StageII_one_against_All[which(colData(vst_Stage_II_sub)$stages=="Stage II")]<-"Stage II"
colData(vst_Stage_II_sub)$StageII_one_against_All[which(colData(vst_Stage_II_sub)$stages!="Stage II")]<-"All others"
colData(vst_Stage_III_sub)$StageIII_one_against_All[which(colData(vst_Stage_III_sub)$stages=="Stage III")]<-"Stage III"
colData(vst_Stage_III_sub)$StageIII_one_against_All[which(colData(vst_Stage_III_sub)$stages!="Stage III")]<-"All others"
#####################################################################################################################
pca_stageI_one_against_All<-plotPCA(vst_Stage_I_sub, intgroup="StageI_one_against_All") + theme_bw() + ggtitle("DE Genes from Stage I")
pca_stageII_one_against_All<-plotPCA(vst_Stage_II_sub, intgroup="StageII_one_against_All") + theme_bw() + ggtitle("DE Genes from Stage II")
pca_stageIII_one_against_All<-plotPCA(vst_Stage_III_sub, intgroup="StageIII_one_against_All") + theme_bw() + ggtitle("DE Genes from Stage III")

pca_plots<-grid.arrange(pca_stageI_one_against_All, pca_stageII_one_against_All,pca_stageIII_one_against_All, nrow = 2)

# FindClusters_resolution
png(filename=paste(output_dir,"Stage_one_against_All.png",sep=""), width = 36, height = 48, res=600, units = "cm")
	pca_plots<-grid.arrange(pca_stageI_one_against_All, pca_stageII_one_against_All,pca_stageIII_one_against_All, nrow = 2)
dev.off()
#####################################################################################################################
Stage_I_sub   <-Stage_I_sub_bck
Stage_II_sub  <-Stage_II_sub_bck
Stage_III_sub <-Stage_III_sub_bck

colData(vst_Stage_I_sub)$Primary_TumorStage_I        <-""
colData(vst_Stage_II_sub)$Primary_TumorStage_II      <-""
colData(vst_Stage_III_sub)$Primary_TumorStage_III    <-""

colData(vst_Stage_I_sub)$Primary_TumorStage_I[which(colData(vst_Stage_I_sub)$Tumor_Stage=="Primary TumorStage I")]<-"Primary Tumor Stage I"
colData(vst_Stage_I_sub)$Primary_TumorStage_I[which(colData(vst_Stage_I_sub)$Tumor_Stage!="Primary TumorStage I")]<-"All others"
colData(vst_Stage_II_sub)$Primary_TumorStage_II[which(colData(vst_Stage_II_sub)$Tumor_Stage=="Primary TumorStage II")]<-"Primary Tumor Stage II"
colData(vst_Stage_II_sub)$Primary_TumorStage_II[which(colData(vst_Stage_II_sub)$Tumor_Stage!="Primary TumorStage II")]<-"All others"
colData(vst_Stage_III_sub)$Primary_TumorStage_III[which(colData(vst_Stage_III_sub)$Tumor_Stage=="Primary TumorStage III")]<-"Primary Tumor Stage III"
colData(vst_Stage_III_sub)$Primary_TumorStage_III[which(colData(vst_Stage_III_sub)$Tumor_Stage!="Primary TumorStage III")]<-"All others"

pca_stageI<-plotPCA(vst_Stage_I_sub, intgroup="Primary_TumorStage_I") + theme_bw() + ggtitle("DE Genes from Primary Tumor_Stage I")
pca_stageII<-plotPCA(vst_Stage_II_sub, intgroup="Primary_TumorStage_II") + theme_bw() + ggtitle("DE Genes from Primary Tumor_Stage II")
pca_stageIII<-plotPCA(vst_Stage_III_sub, intgroup="Primary_TumorStage_III") + theme_bw() + ggtitle("DE Genes from Primary Tumor_Stage III")

# FindClusters_resolution
png(filename=paste(output_dir,"Stage_PrimaryTumorStage_against_All.png",sep=""), width = 36, height = 48, res=600, units = "cm")
	plot_grid(pca_stageI, pca_stageII,pca_stageIII, ncol = 2, nrow = 3)
dev.off()
#####################################################################################################################
Stage_I_sub   <-Stage_I_sub_bck
Stage_II_sub  <-Stage_II_sub_bck
Stage_III_sub <-Stage_III_sub_bck

# Primary Tumor Stage I Normal
# Primary Tumor Stage I Tumor
# All others : colour black
colData(vst_Stage_I_sub)$Type_Stage_Tumor_Stage_I    <-""
colData(vst_Stage_I_sub)$Type_Stage_Tumor_Stage_II   <-""
colData(vst_Stage_I_sub)$Type_Stage_Tumor_Stage_III   <-""

which(colData(vst_Stage_I_sub)$Tumor_Stage!="Primary TumorStage I")

########################################################################################################################
library(tidyverse)
library(ggrepel)
library(kableExtra)

# Sort table by abs(log2FoldChange) and -log(padj)
Normal_Tumor_sort<- Normal_Tumor[order(Normal_Tumor$padj), ]

# Field for top 2.5 percent of sorted sample
Normal_Tumor_sort$Normal_Tumor_sort_2.5<-FALSE

# Field for top 2.5 percent of sorted sample
Normal_Tumor_sort[1:(dim(Normal_Tumor_sort)[1]*0.025),"Normal_Tumor_sort_2.5"]<-TRUE

# "Unchanged"
Normal_Tumor_sort$Expression<-0

# Set expression
Normal_Tumor_sort[intersect(which(Normal_Tumor_sort$Normal_Tumor_sort_2.5), which(Normal_Tumor_sort$log2FoldChange < 0)),"Expression"]<--1
Normal_Tumor_sort[intersect(which(Normal_Tumor_sort$Normal_Tumor_sort_2.5), which(Normal_Tumor_sort$log2FoldChange >= 0)),"Expression"]<-1

Normal_Tumor_sort$Categories<-""
Normal_Tumor_sort[which(Normal_Tumor_sort$Expression==0),"Categories"]<-"Uncategorized"
Normal_Tumor_sort[which(Normal_Tumor_sort$Expression==1),"Categories"]<-"Up-regulated"
Normal_Tumor_sort[which(Normal_Tumor_sort$Expression==-1),"Categories"]<-"Down-regulated"

# Create volcano plot
p1 <- ggplot(Normal_Tumor_sort, aes(log2FoldChange, -log(padj))) + # -log10 conversion  
  geom_point(size = 2/5) +  theme_bw()

# Adding color to differentially expressed genes (DEGs)
p2 <- ggplot(Normal_Tumor_sort, aes(log2FoldChange, -log(padj),color = Categories)) + geom_point(size = 2/5,aes(color = Categories))  +
  xlab(expression("log2FoldChange")) + 
  ylab(expression("-log"[10]*"padj")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_bw() + ggtitle(paste("2.5% DE Genes Tumor vs. Normal Samples \nsorted by padj\n",paste("Up-regulated :",sum(Normal_Tumor_sort$Categories=="Up-regulated"),"Down-regulated :",sum(Normal_Tumor_sort$Categories=="Down-regulated"),sep=" ")))
########################################################################################################################
# FindClusters_resolution
png(filename=paste(output_dir,"Volcano_Plot_Normal_Tumor.png",sep=""), width = 16, height = 16, res=600, units = "cm")
	p2
dev.off()



			 
