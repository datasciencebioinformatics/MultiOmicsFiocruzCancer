# Reat stage specific genes 
stage_specific_genes <-  read.xlsx(file="/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Table3.xlsx", 2)   # read first sheet

stage_specific_genes$stage<-""
stage_specific_genes[rownames(stage_specific_genes) %in% unique_stage_I,"stage"]<-"Stage I"
stage_specific_genes[rownames(stage_specific_genes) %in% unique_stage_II,"stage"]<-"Stage II"
stage_specific_genes[rownames(stage_specific_genes) %in% unique_stage_III,"stage"]<-"Stage III"

write_tsv(df_mean, paste(output_dir,"Table2.tsv",sep=""))			




# Genes that are tumor genes (RPKM ≥ 4, tumor vs. normal samples : log2foldchange >1 and paired t-test FDR ≤ 0.05) but also stage-specific genes (RPKM ≥ 4, stage-specific vs. normal samples : log2foldchange >1 and paired t-test FDR ≤ 0.05).
genes_stages_I   <-list_per_stage_comparisson$sample_stage_I_sample_stage_II[which(list_per_stage_comparisson$sample_stage_I_sample_stage_II$log2change>0.0 & list_per_stage_comparisson$sample_stage_I_sample_stage_III$log2change>0.0),"gene"]
genes_stages_II  <-list_per_stage_comparisson$sample_stage_II_sample_stage_I[which(list_per_stage_comparisson$sample_stage_II_sample_stage_I$log2change>0.20 & list_per_stage_comparisson$sample_stage_II_sample_stage_III$log2change>0.20),"gene"]
genes_stages_III <-list_per_stage_comparisson$sample_stage_III_sample_stage_I[which(list_per_stage_comparisson$sample_stage_III_sample_stage_I$log2change>0.30 & list_per_stage_comparisson$sample_stage_III_sample_stage_II$log2change>0.30),"gene"]

genes_stages_I<-genes_stages_I[genes_stages_I %in% genes_Stage_I$gene]
genes_stages_II<-genes_stages_II[genes_stages_II %in% genes_Stage_II$gene]
genes_stages_III<-genes_stages_III[genes_stages_III %in% genes_Stage_III$gene]

p_stage_I_unpaired<-ggplot(unstranded_data_samples_unapaired[unstranded_data_samples_unapaired$gene_id %in% genes_stages_I,], aes(x=stages, y=RPKM, fill=stages)) +  geom_boxplot()+ facet_wrap(~SYMBOL, nrow = 5,ncol = 6, scales="free")+ theme_bw() + ggtitle("Biomarkers for stage I") + scale_fill_manual(values=c('#e8f2a1', '#729fcf', '#ffaa95')) + stat_summary(fun.y=mean, fun.min = min, fun.max = max ,geom="point", shape=20, size=8, color="red", fill="red")+ geom_jitter(aes(colour = stages), size=0.4, alpha=0.9) + scale_color_manual(values=c('#e8f2a1', '#729fcf', '#ffaa95'))
p_stage_II_unpaired<-ggplot(unstranded_data_samples_unapaired[unstranded_data_samples_unapaired$gene_id %in% genes_stages_II,], aes(x=stages, y=RPKM, fill=stages)) +  geom_boxplot()+ facet_wrap(~SYMBOL, nrow = 5,ncol = 6, scales="free")+ theme_bw() + ggtitle("Biomarkers for stage II") + scale_fill_manual(values=c('#e8f2a1', '#729fcf', '#ffaa95')) + stat_summary(fun.y=mean, fun.min = min, fun.max = max ,geom="point", shape=20, size=8, color="red", fill="red")+ geom_jitter(aes(colour = stages), size=0.4, alpha=0.9) + scale_color_manual(values=c('#e8f2a1', '#729fcf', '#ffaa95'))
p_stage_III_unpaired<-ggplot(unstranded_data_samples_unapaired[unstranded_data_samples_unapaired$gene_id %in% genes_stages_III,], aes(x=stages, y=RPKM, fill=stages)) +  geom_boxplot()+ facet_wrap(~SYMBOL, nrow = 5,ncol = 6, scales="free")+ theme_bw() + ggtitle("Biomarkers for stage III") + scale_fill_manual(values=c('#e8f2a1', '#729fcf', '#ffaa95')) + stat_summary(fun.y=mean, fun.min = min, fun.max = max ,geom="point", shape=20, size=8, color="red", fill="red")+ geom_jitter(aes(colour = stages), size=0.4, alpha=0.9) + scale_color_manual(values=c('#e8f2a1', '#729fcf', '#ffaa95'))

# p_stage_I_unpaired.png
png(filename=paste(output_folder,paste("p_stage_I_unpaired.png",sep=""),sep=""), width = 24, height = 24, res=600, units = "cm")
  print(p_stage_I_unpaired + theme(legend.position="bottom"))
dev.off()

# p_stage_II_unpaired.png
png(filename=paste(output_folder,paste("p_stage_II_unpaired.png",sep=""),sep=""), width = 24, height = 24, res=600, units = "cm")
  print(p_stage_II_unpaired + theme(legend.position="bottom"))
dev.off()

# p_stage_III_unpaired.png
png(filename=paste(output_folder,paste("p_stage_III_unpaired.png",sep=""),sep=""), width = 24, height = 24, res=600, units = "cm")
  print(p_stage_III_unpaired + theme(legend.position="bottom"))
dev.off()


biomarkers<-c("ABCA8" , "ADAMTS8" , "ALK" , "APC" , "ARHGEF12" , "ASPA" , "AURKA" , "AURKB" , "BIRC5" , "BRAF" , "BRCA" , "CBFA2T3" , "CCNB1" , "CCNB2" , "CDK1" , "CDKN1C" , "CEP55" , "CHEK1" , "CKIT" , "CMET " , "DAB2IP" , "DCC" , "DDR2" , "DDX5" , "EGFR" , "ERBB2" , "ERCC1" , "EXT1" , "FGFR" , "FHL1" , "FOXP1" , "GPC3" , "HER2" , "JAK2 " , "KCNRG" , "KLK10" , "KRAS" , "LATS2" , "LIMD1" , "MCC" , "MET" , "NBL1" , "NTRK1" , "P53" , "PCNA" , "PD-1 " , "PD-L1" , "PD-L2 " , "PIK3CA" , "PTCH1" , "PTEN" , "PYCR1" , "RAMP3" , "RAP1A" , "RASSF2" , "RB" , "RECK" , "RET" , "RHOB" , "ROS" , "RRM1" , "SASH1" , "STARD13" , "TBRG1" , "TKIS " , "TOP2A" , "TP53" , "TPX2" , "TRIM13" , "TYMS" , "UBE2C" )
p_stage_biomarkers_unpaired<-ggplot(unstranded_data_samples_unapaired[unstranded_data_samples_unapaired$SYMBOL %in% biomarkers,], aes(x=stages, y=RPKM, fill=stages)) +  geom_boxplot()+ facet_wrap(~SYMBOL, nrow = 5,ncol = 6, scales="free")+ theme_bw() + scale_fill_manual(values=c('#e8f2a1', '#729fcf', '#ffaa95')) + stat_summary(fun.y=mean, fun.min = min, fun.max = max ,geom="point", shape=20, size=8, color="red", fill="red")+ geom_jitter(aes(colour = stages), size=0.4, alpha=0.9) + scale_color_manual(values=c('#e8f2a1', '#729fcf', '#ffaa95')) + theme(axis.text.x = element_text(angle = 90)) 
p_stage_biomarkers_paired<-ggplot(unstranded_data_samples[unstranded_data_samples$SYMBOL %in% biomarkers,], aes(x=tissue_type, y=RPKM, fill=tissue_type)) +  geom_boxplot()+ facet_wrap(~SYMBOL, nrow = 5,ncol = 6, scales="free")+ theme_bw() + scale_fill_manual(values=c('#e8f2a1', '#729fcf', '#ffaa95')) + stat_summary(fun.y=mean, fun.min = min, fun.max = max ,geom="point", shape=20, size=8, color="red", fill="red")+ geom_jitter(aes(colour = stages), size=0.4, alpha=0.9) + scale_color_manual(values=c('#e8f2a1', '#729fcf', '#ffaa95'))+ theme(axis.text.x = element_text(angle = 90)) 

# p_stage_III_unpaired.png
png(filename=paste(output_folder,paste("p_stage_biomarkers_unpaired.png",sep=""),sep=""), width = 24, height = 24, res=600, units = "cm")
  print(p_stage_biomarkers_unpaired + theme(legend.position="bottom"))
dev.off()

# p_stage_III_unpaired.png
png(filename=paste(output_folder,paste("p_stage_biomarkers_paired.png",sep=""),sep=""), width = 24, height = 24, res=600, units = "cm")
  print(p_stage_biomarkers_paired + theme(legend.position="bottom"))
dev.off()







#######################################################################################################################
biomarkers<-c("AMTN" , "FABP7" , "OLFM4")
# Complete dataset
unstranded_data_samples_complete<-rbind(unstranded_data_samples_unapaired[,c("SYMBOL","RPKM","stages","tissue_type")],unstranded_data_samples[unstranded_data_samples$tissue_type=="Normal",c("SYMBOL","RPKM","stages","tissue_type")])

biomarkers_AMTN_boxplot2   <-ggplot(unstranded_data_samples_complete[unstranded_data_samples_complete$SYMBOL %in% biomarkers,], aes(x=stages, y=RPKM, fill=tissue_type, color=tissue_type)) + facet_wrap(~SYMBOL, nrow = 5,ncol = 6, scales="free")+ theme_bw() + theme(axis.text.x = element_text(angle = 90)) +  geom_jitter(width = 0.25) 
biomarkers_AMTN_errobar2   <-ggplot(unstranded_data_samples_complete[unstranded_data_samples_complete$SYMBOL %in% biomarkers,], aes(x=stages, y=RPKM, fill=tissue_type, color=tissue_type)) + facet_wrap(~SYMBOL, nrow = 5,ncol = 6, scales="free")+ theme_bw() + theme(axis.text.x = element_text(angle = 90)) +  stat_boxplot(geom ='errorbar', width = 0.5) 



# p_stage_III_unpaired.png
png(filename=paste(output_folder,paste("biomarkers_AMTN_boxplot2.png",sep=""),sep=""), width = 16, height = 8, res=600, units = "cm")
  print(biomarkers_AMTN_boxplot2 + theme(legend.position="bottom"))
dev.off()

# p_stage_III_unpaired.png
png(filename=paste(output_folder,paste("biomarkers_AMTN_errobar2.png",sep=""),sep=""), width = 16, height = 8, res=600, units = "cm")
  print(biomarkers_AMTN_errobar2 + theme(legend.position="bottom"))
dev.off()



#######################################################################################################################






#######################################################################################################################
# Reload colData from file
# Reload unstranded_data from file
###########################################################################################################################
merged_data_patient_info_file       <- "/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv"                  #
colData_file                        <- "/home/felipe/Documentos/LungPortal/samples/colData.tsv"                           #
###########################################################################################################################
unstranded_data                    <-unstranded_data_filter
merged_data_patient_info_data      <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)#
colData_data                       <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)                 #
rownames(colData)                  <-colData$patient_id                                                                   #
###########################################################################################################################
#omit NA values from vector
unstranded_data <- na.omit(unstranded_data)
########################################################################################################################
# A panel to analyse differential Category comparing samples of each stage against all others stages.
########################################################################################################################
# Only tumor samples
colData_tumor <-colData[colData$tissue_type=="Tumor",]
colData_normal<-colData[colData$tissue_type=="Normal",]

# Vector with each stage
stages_str<-c("stage_I","stage_II","stage_III")

# Samples of each stage stored in colData                                                                                             #
sample_stage_I  <-colData_tumor[colData_tumor$stages=="Stage I","patient_id"]                                                                     #
sample_stage_II <-colData_tumor[colData_tumor$stages=="Stage II","patient_id"]                                                                    #
sample_stage_III<-colData_tumor[colData_tumor$stages=="Stage III","patient_id"]                                                                   #
sample_normal   <-colData_normal[,"patient_id"]                                                                   #
#######################################################################################################################################
df_table_comparisson=rbind(data.frame(Stage_i="sample_stage_I",Stage_ii="sample_normal"),
data.frame(Stage_i="sample_stage_II",Stage_ii="sample_normal"),
data.frame(Stage_i="sample_stage_III",Stage_ii="sample_normal"))
####################################################################################################################
# Take p-value
df_mean<-data.frame(ENSEMBL=stage_specific_genes$ENSEMBL,ENTREZID=stage_specific_genes$ENTREZID,SYMBOL=stage_specific_genes$SYMBOL,
avg.normal=rowMeans(unstranded_data[stage_specific_genes$gene,sample_normal]),
std.normal=0,avg.stageI=rowMeans(unstranded_data[stage_specific_genes$gene,sample_stage_I]),
std.stageI=0, avg.stageII=rowMeans(unstranded_data[stage_specific_genes$gene,sample_stage_II]), std.stageII=0, avg.stageIII=rowMeans(unstranded_data[stage_specific_genes$gene,sample_stage_III]), std.stageIII=0)

# For each gene, calculate too the 
for (gene in rownames(df_mean))
{
  df_mean[gene,"std.normal"]<-sd(unstranded_data[gene,sample_normal])
  df_mean[gene,"std.stageI"]<-sd(unstranded_data[gene,sample_stage_I])
  df_mean[gene,"std.stageII"]<-sd(unstranded_data[gene,sample_stage_II])
  df_mean[gene,"std.stageIII"]<-sd(unstranded_data[gene,sample_stage_III])  
}
df_mean$stage<-""
df_mean[rownames(df_mean) %in% unique_stage_I,"stage"]<-"Stage I"
df_mean[rownames(df_mean) %in% unique_stage_II,"stage"]<-"Stage II"
df_mean[rownames(df_mean) %in% unique_stage_III,"stage"]<-"Stage III"

write_tsv(df_mean, paste(output_dir,"Table3.tsv",sep=""))			
