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

biomarkers_AMTN_boxplot2   <-ggplot(unstranded_data_samples_complete[unstranded_data_samples_complete$SYMBOL %in% selected_genes,], aes(x=stages, y=RPKM, fill=tissue_type, color=tissue_type)) + facet_wrap(~SYMBOL, nrow = 5,ncol = 6, scales="free")+ theme_bw() + theme(axis.text.x = element_text(angle = 90)) +  geom_jitter(width = 0.25) 
biomarkers_AMTN_errobar2   <-ggplot(unstranded_data_samples_complete[unstranded_data_samples_complete$SYMBOL %in% biomarkers,], aes(x=stages, y=RPKM, fill=tissue_type, color=tissue_type)) + facet_wrap(~SYMBOL, nrow = 5,ncol = 6, scales="free")+ theme_bw() + theme(axis.text.x = element_text(angle = 90)) +  stat_boxplot(geom ='errorbar', width = 0.5) 



# p_stage_III_unpaired.png
png(filename=paste(output_folder,paste("biomarkers_AMTN_boxplot2.png",sep=""),sep=""), width = 24, height = 10, res=500, units = "cm")
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










####################################################################################################################
# Take p-value
selected_genes_Stage_I_data$avg.normal<-rowMeans(unstranded_data[selected_genes_Stage_I_data$gene,sample_normal])
selected_genes_Stage_I_data$std.normal<-0
selected_genes_Stage_I_data$avg.stageI<-rowMeans(unstranded_data[selected_genes_Stage_I_data$gene,sample_stage_I])
selected_genes_Stage_I_data$std.stageI<-0

selected_genes_Stage_II_data$avg.normal<-rowMeans(unstranded_data[selected_genes_Stage_II_data$gene,sample_normal])
selected_genes_Stage_II_data$std.normal<-0
selected_genes_Stage_II_data$avg.stageII<-rowMeans(unstranded_data[selected_genes_Stage_II_data$gene,sample_stage_II])
selected_genes_Stage_II_data$std.stageII<-0

selected_genes_Stage_III_data$avg.normal<-rowMeans(unstranded_data[selected_genes_Stage_III_data$gene,sample_normal])
selected_genes_Stage_III_data$std.normal<-0
selected_genes_Stage_III_data$avg.stageIII<-rowMeans(unstranded_data[selected_genes_Stage_III_data$gene,sample_stage_III])
selected_genes_Stage_III_data$std.stageIII<-0


# For each gene, calculate too the 
for (gene in rownames(selected_genes_Stage_I_data))
{
  selected_genes_Stage_I_data[gene,"std.normal"]<-sd(unstranded_data[gene,sample_normal])
  selected_genes_Stage_I_data[gene,"std.stageI"]<-sd(unstranded_data[gene,sample_stage_I])  
}

# For each gene, calculate too the 
for (gene in rownames(selected_genes_Stage_II_data))
{
  selected_genes_Stage_II_data[gene,"std.normal"]<-sd(unstranded_data[gene,sample_normal])
  selected_genes_Stage_II_data[gene,"std.stageII"]<-sd(unstranded_data[gene,sample_stage_II])  
}

# For each gene, calculate too the 
for (gene in rownames(selected_genes_Stage_III_data))
{
  selected_genes_Stage_III_data[gene,"std.normal"]<-sd(unstranded_data[gene,sample_normal])
  selected_genes_Stage_III_data[gene,"std.stageIII"]<-sd(unstranded_data[gene,sample_stage_III])  
}

#############
# For each gene, add gene_id
selected_genes_Stage_I_data$gene_id<-""
for (gene_row in rownames(selected_genes_Stage_I_data))
{	
	# Store gene id in the vector
	# Simply trim the gene id before the "." to save it in the ENSEML format
	selected_genes_Stage_I_data[gene_row,"gene_id"]<-strsplit(selected_genes_Stage_I_data[gene_row,"gene"], split = "\\.")[[1]][1]	
}

# For each gene, add gene_id
selected_genes_Stage_II_data$gene_id<-""
for (gene_row in rownames(selected_genes_Stage_II_data))
{	
	# Store gene id in the vector
	# Simply trim the gene id before the "." to save it in the ENSEML format
	selected_genes_Stage_II_data[gene_row,"gene_id"]<-strsplit(selected_genes_Stage_II_data[gene_row,"gene"], split = "\\.")[[1]][1]	
}

# For each gene, add gene_id
selected_genes_Stage_III_data$gene_id<-""
for (gene_row in rownames(selected_genes_Stage_III_data))
{	
	# Store gene id in the vector
	# Simply trim the gene id before the "." to save it in the ENSEML format
	selected_genes_Stage_III_data[gene_row,"gene_id"]<-strsplit(selected_genes_Stage_III_data[gene_row,"gene"], split = "\\.")[[1]][1]	
}

# ids_stage_I - all ENSEMBL anotated using bitr
ids_stage_I       <-bitr(selected_genes_Stage_I_data$gene_id, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
ids_stage_II      <-bitr(selected_genes_Stage_II_data$gene_id, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
ids_stage_III     <-bitr(selected_genes_Stage_III_data$gene_id, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")

colnames(selected_genes_Stage_I_data)[10]<-"ENSEMBL"
colnames(selected_genes_Stage_II_data)[10]<-"ENSEMBL"
colnames(selected_genes_Stage_III_data)[10]<-"ENSEMBL"

selected_genes_Stage_I_data<-merge(selected_genes_Stage_I_data,ids_stage_I,by="ENSEMBL")
selected_genes_Stage_II_data<-merge(selected_genes_Stage_II_data,ids_stage_II,by="ENSEMBL")
selected_genes_Stage_III_data<-merge(selected_genes_Stage_III_data,ids_stage_III,by="ENSEMBL")

write_tsv(selected_genes_Stage_I_data, paste(output_dir,"all_genes_Stage_I_data.tsv",sep=""))			
write_tsv(selected_genes_Stage_II_data, paste(output_dir,"all_genes_Stage_II_data.tsv",sep=""))	
write_tsv(selected_genes_Stage_III_data, paste(output_dir,"all_genes_Stage_III_data.tsv",sep=""))		

selected_genes_Stage_I_data_bck<-selected_genes_Stage_I_data[selected_genes_Stage_I_data$avg.normal<2,]
selected_genes_Stage_I_data_bck$lfc<-(selected_genes_Stage_I_data_bck$avg.stageI/selected_genes_Stage_I_data_bck$avg.normal)
selected_genes_Stage_I_data_bck<- selected_genes_Stage_I_data_bck[selected_genes_Stage_I_data_bck$lfc>100,]

selected_genes_Stage_II_data_bck<-selected_genes_Stage_II_data[selected_genes_Stage_II_data$avg.normal<2,]
selected_genes_Stage_II_data_bck$lfc<-(selected_genes_Stage_II_data_bck$avg.stageII/selected_genes_Stage_II_data_bck$avg.normal)
selected_genes_Stage_II_data_bck<- selected_genes_Stage_II_data_bck[selected_genes_Stage_II_data_bck$lfc>100,]

selected_genes_Stage_III_data_bck<-selected_genes_Stage_III_data[selected_genes_Stage_III_data$avg.normal<2,]
selected_genes_Stage_III_data_bck$lfc<-(selected_genes_Stage_III_data_bck$avg.stageIII/selected_genes_Stage_III_data_bck$avg.normal)
selected_genes_Stage_III_data_bck<- selected_genes_Stage_III_data_bck[selected_genes_Stage_III_data_bck$lfc>100,]

unique_stage_I  =intersect(setdiff(selected_genes_Stage_I_data_bck$SYMBOL, c(selected_genes_Stage_II_data_bck$SYMBOL,selected_genes_Stage_III_data_bck$SYMBOL)),selected_genes_Stage_I_data_bck$SYMBOL)
unique_stage_II  =intersect(setdiff(selected_genes_Stage_II_data_bck$SYMBOL, c(selected_genes_Stage_I_data_bck$SYMBOL,selected_genes_Stage_III_data_bck$SYMBOL)),selected_genes_Stage_II_data_bck$SYMBOL)
unique_stage_III  =intersect(setdiff(selected_genes_Stage_III_data_bck$SYMBOL, c(selected_genes_Stage_I_data_bck$SYMBOL,selected_genes_Stage_II_data_bck$SYMBOL)),selected_genes_Stage_III_data_bck$SYMBOL)

selected_genes_Stage_I_data_bck <-selected_genes_Stage_I_data_bck[selected_genes_Stage_I_data_bck$SYMBOL %in% unique_stage_I,]
selected_genes_Stage_II_data_bck<-selected_genes_Stage_II_data_bck[selected_genes_Stage_II_data_bck$SYMBOL %in% unique_stage_II,]
selected_genes_Stage_III_data_bck<-selected_genes_Stage_III_data_bck[selected_genes_Stage_III_data_bck$SYMBOL %in% unique_stage_III,]
	
df_stage_I<-data.frame(SYMBOL=selected_genes_Stage_I_data_bck$SYMBOL,
avg.normal=selected_genes_Stage_I_data_bck$avg.normal,
std.normal=selected_genes_Stage_I_data_bck$std.normal,
avg.stageI=selected_genes_Stage_I_data_bck$avg.stageI,
std.stageI=selected_genes_Stage_I_data_bck$std.stageI,
FC=selected_genes_Stage_I_data_bck$lfc,stage="Stage I")

df_stage_II<-data.frame(SYMBOL=selected_genes_Stage_II_data_bck$SYMBOL,
avg.normal=selected_genes_Stage_II_data_bck$avg.normal,
std.normal=selected_genes_Stage_II_data_bck$std.normal,
avg.stageI=selected_genes_Stage_II_data_bck$avg.stageII,
std.stageI=selected_genes_Stage_II_data_bck$std.stageII,
FC=selected_genes_Stage_II_data_bck$lfc,stage="Stage II")

df_stage_III<-data.frame(SYMBOL=selected_genes_Stage_III_data_bck$SYMBOL,
avg.normal=selected_genes_Stage_III_data_bck$avg.normal,
std.normal=selected_genes_Stage_III_data_bck$std.normal,
avg.stageI=selected_genes_Stage_III_data_bck$avg.stageIII,
std.stageI=selected_genes_Stage_III_data_bck$std.stageIII,
FC=selected_genes_Stage_III_data_bck$lfc,stage="Stage III")

rbind(df_stage_I,df_stage_II,df_stage_III)


####################################################################################################################
table_I<-c("KRT14", "KRT16", "NTS", "SPRR1B", "GPX2", "SPRR1B", "AKR1B10", "GPX2", "AKR1B10", "KRT13", "SPRR2A", "KRT13", "AKR1B10", "KRT6B", "S100A7")
table_II<-c("GRB7", "SRCv", "RNPS1v", "HOOK2", "EFTUD2","PRKCI","DVL2","HAUS1","RNF2","PHB1","ELOC","PSMC6","THAP7","SEH1L")
table_III<-c("CDK8","KRT1","NEDD1","GMCL1","GOLT1B","BEX2","PRMT6","RBBP7","SCNM1","TP53","CEP131","CLK2","EHMT2","FOXK2","PNKP","PRMT5","USP21")
table_IV<-c("MAGEA6","KRT31","KRT75","KRT16","FOXE1","CRCT1","PITX1","KRT15","TP63","TFAP2A","NUF2","FOXM1","ANLN","BUB1B","CEP55","PLK1")

selected_genes<-unique(c(table_I,table_II,table_III,table_IV))

selected_genes<-unique(c(table_I))


# Control samples df
df_control_samples<-unstranded_data_samples[unstranded_data_samples$tissue_type=="Normal",c("SYMBOL","RPKM","stages","tissue_type")]
df_control_samples$stages<-"Control"
df_control_samples$stages= factor(df_control_samples$stages, levels=c("Stage I","Stage II","Stage III","Control"))

unstranded_data_samples_complete<-rbind(unstranded_data_samples_unapaired[,c("SYMBOL","RPKM","stages","tissue_type")],df_control_samples)
biomarkers_AMTN_boxplot2   <-ggplot(unstranded_data_samples_complete[unstranded_data_samples_complete$SYMBOL %in% selected_genes,], aes(x=stages, y=RPKM, fill=stages, color=stages)) + facet_wrap(~SYMBOL, nrow = 3,ncol = 5, scales="free")+ theme_bw() + theme(axis.text.x = element_text(angle = 90)) +  geom_boxplot(outlier.shape=NA) +  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_fill_manual(values=c('black','#e8f2a1', '#729fcf', '#ffaa95')) + scale_colour_manual(values=c('black','#e8f2a1', '#729fcf', '#ffaa95')) + theme(legend.position="bottom")+ geom_jitter(aes(colour = stages), size=0.4, alpha=0.9)


# p_stage_III_unpaired.png
png(filename=paste(output_folder,paste("Figure2.png",sep=""),sep=""), width = 16, height = 10, res=600, units = "cm")
  print(biomarkers_AMTN_boxplot2 + theme(legend.position="bottom"))
dev.off()







####################################################################################################################
stage_I_GO_GENES<-unique(df_all_annotation[df_all_annotation$CluterProfiler %in% stage_I_GO,"Symbol"])
stage_II_GO_GENES<-unique(df_all_annotation[df_all_annotation$CluterProfiler %in% stage_II_GO,"Symbol"])
stage_III_GO_GENES<-unique(df_all_annotation[df_all_annotation$CluterProfiler %in% stage_III_GO,"Symbol"])

stage_I_Reactome_GENES<-unique(df_all_annotation[df_all_annotation$CluterProfiler %in% stage_I_Reactome,"Symbol"])
stage_II_Reactome_GENES<-unique(df_all_annotation[df_all_annotation$CluterProfiler %in% stage_II_Reactome,"Symbol"])
stage_III_Reactome_GENES<-unique(df_all_annotation[df_all_annotation$CluterProfiler %in% stage_III_Reactome,"Symbol"])

stage_I_KEGG_GENES<-unique(df_all_annotation[df_all_annotation$CluterProfiler %in% stage_I_KEGG,"Symbol"])
stage_II_KEGG_Genes<-unique(df_all_annotation[df_all_annotation$CluterProfiler %in% stage_II_KEGG,"Symbol"])
stage_III_KEGG_Genes<-unique(df_all_annotation[df_all_annotation$CluterProfiler %in% stage_III_KEGG,"Symbol"])
####################################################################################################################
table_I[table_I %in% c(stage_I_GO_GENES,stage_I_Reactome_GENES,stage_I_KEGG_GENES)]
table_I[table_I %in% c(stage_II_GO_GENES,stage_II_Reactome_GENES,stage_II_KEGG_Genes)]
table_I[table_I %in% c(stage_III_GO_GENES,stage_III_Reactome_GENES,stage_III_KEGG_Genes)]

table_II[table_II %in% c(stage_I_GO_GENES,stage_I_Reactome_GENES,stage_I_KEGG_GENES)]
table_II[table_II %in% c(stage_II_GO_GENES,stage_II_Reactome_GENES,stage_II_KEGG_Genes)]
table_II[table_II %in% c(stage_III_GO_GENES,stage_III_Reactome_GENES,stage_III_KEGG_Genes)]

table_III[table_III %in% c(stage_I_GO_GENES,stage_I_Reactome_GENES,stage_I_KEGG_GENES)]
table_III[table_III %in% c(stage_II_GO_GENES,stage_II_Reactome_GENES,stage_II_KEGG_Genes)]
table_III[table_III %in% c(stage_III_GO_GENES,stage_III_Reactome_GENES,stage_III_KEGG_Genes)]

table_IV[table_IV %in% c(stage_I_GO_GENES,stage_I_Reactome_GENES,stage_I_KEGG_GENES)]
table_IV[table_IV %in% c(stage_II_GO_GENES,stage_II_Reactome_GENES,stage_II_KEGG_Genes)]
table_IV[table_IV %in% c(stage_III_GO_GENES,stage_III_Reactome_GENES,stage_III_KEGG_Genes)]













####################################################################################################################
biomarkers<-c("PDCD11" )
unstranded_data_samples_complete<-rbind(unstranded_data_samples_unapaired[,c("SYMBOL","RPKM","stages","tissue_type")],df_control_samples)
biomarkers_AMTN_boxplot2   <-ggplot(unstranded_data_samples_complete[unstranded_data_samples_complete$SYMBOL %in% biomarkers,], aes(x=stages, y=RPKM, fill=stages, color=stages)) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +  geom_boxplot(outlier.shape=NA) +  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_fill_manual(values=c('black','#e8f2a1', '#729fcf', '#ffaa95')) + scale_colour_manual(values=c('black','#e8f2a1', '#729fcf', '#ffaa95')) + theme(legend.position="bottom")+ geom_jitter(aes(colour = stages), size=0.4, alpha=0.9) + ggtile("PDCD11")

# p_stage_III_unpaired.png
png(filename=paste(output_folder,paste("p_stage_biomarkers_pdcd11.png",sep=""),sep=""), width = 12, height = 12, res=600, units = "cm")
  print(biomarkers_AMTN_boxplot2 + theme(legend.position="bottom"))
dev.off()

# p_stage_III_unpaired.png
png(filename=paste(output_folder,paste("p_stage_biomarkers_paired.png",sep=""),sep=""), width = 24, height = 24, res=600, units = "cm")
  print(p_stage_biomarkers_paired + theme(legend.position="bottom"))
dev.off()



####################################################################################################################
biomarkers<-c("ENSG00000148843")

