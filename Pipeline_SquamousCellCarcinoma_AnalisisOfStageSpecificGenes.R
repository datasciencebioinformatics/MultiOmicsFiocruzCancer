# Reat stage specific genes 
stage_specific_genes <-  read.xlsx(file="/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Table3.xlsx", 2)   # read first sheet

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


stage_specific_genes
