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
rownames_file_2         <- "/home/felipe/Documentos/LungPortal/samples/gene_ids.txt"

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
/home/felipe/Documentos/LungPortal/Table_1.xls
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

# If you use a design of ~0 + covariable, then you will have a coefficient for each level of condition. 
# This is one of the cases where it helps to not have condition at the end of the design (for convenience, we often recommend to put condition at the end, but not in this case).

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
pca_tumor_stages<-plotPCA(vst_dds, intgroup="Tumor_Stage")+ theme_bw()         + ggtitle("stages")

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
Stage_I    <-data.frame(results(dds,contrast=list(c("stagesStage.I"), c("stagesStage.II","stagesStage.III"))))
Stage_II   <-data.frame(results(dds,contrast=list(c("stagesStage.II"), c("stagesStage.I","stagesStage.III")))) 
Stage_III   <-data.frame(results(dds,contrast=list(c("stagesStage.III"), c("stagesStage.I","stagesStage.II")))) 

log2fc_threshold  = 1
padj_treshold     = 0.05

# Filterby lfcThreshold values
Stage_I_sub<-Stage_I[which(Stage_I$log2FoldChange>log2fc_threshold),]
Stage_II_sub<-Stage_II[which(Stage_II$log2FoldChange>log2fc_threshold),]
Stage_III_sub<-Stage_III[which(Stage_III$log2FoldChange>log2fc_threshold),]

# Save differential expression table
write_tsv(Stage_I_sub, "/home/felipe/Documentos/LungPortal/samples/stages_StageI.tsv")
write_tsv(Stage_II, "/home/felipe/Documentos/LungPortal/samples/stages_StageII.tsv")
write_tsv(Stage_II_sub, "/home/felipe/Documentos/LungPortal/samples/stages_StageIII.tsv")
#####################################################################################################################
