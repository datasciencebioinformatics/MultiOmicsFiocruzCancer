# Reat stage specific genes 
stage_specific_genes <-  read.xlsx(file="/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Table3.xlsx", 2)   # read first sheet

# Genes that are tumor genes (RPKM ≥ 4, tumor vs. normal samples : log2foldchange >1 and paired t-test FDR ≤ 0.05) but also stage-specific genes (RPKM ≥ 4, stage-specific vs. normal samples : log2foldchange >1 and paired t-test FDR ≤ 0.05).

# (RPKM ≥ 4, tumor vs. normal samples : log2foldchange >1 and paired t-test FDR ≤ 0.05)
tumor_genes<-stage_specific_genes[stage_specific_genes$tumor.normal.log2fc >1 & stage_specific_genes$tumor.normal.FDR < 0.05,]

stage_specific_genes_stage_I<-tumor_genes[tumor_genes$gene %in% stage_I_selection$gene,]
stage_specific_genes_stage_II<-tumor_genes[tumor_genes$gene %in% stage_II_selection$gene,]
stage_specific_genes_stage_III<-tumor_genes[tumor_genes$gene %in% stage_III_selection$gene,]

# Stage I, comparisson statistics against Stage II and Stage III
# Stage I   :  Gene is selected if Stage I > Stage II and  Stage I > Stage III
# Stage II  :  Gene is selected if Stage II > Stage I and  Stage II > Stage III
# Stage III : selected if Stage III > Stage I and  Stage III > Stage II
list_per_stage_comparisson$sample_stage_I_sample_stage_II

genes_stages_I   <-list_per_stage_comparisson$sample_stage_I_sample_stage_II[which(list_per_stage_comparisson$sample_stage_I_sample_stage_II$log2change>1 | list_per_stage_comparisson$sample_stage_I_sample_stage_III$log2change>1),"gene"]
genes_stages_II  <-list_per_stage_comparisson$sample_stage_II_sample_stage_I[which(list_per_stage_comparisson$sample_stage_II_sample_stage_I$log2change>1 | list_per_stage_comparisson$sample_stage_II_sample_stage_III$log2change>1),"gene"]
genes_stages_III <-list_per_stage_comparisson$sample_stage_III_sample_stage_I[which(list_per_stage_comparisson$sample_stage_III_sample_stage_I$log2change>1 | list_per_stage_comparisson$sample_stage_III_sample_stage_II$log2change>1),"gene"]

genes_stages_I<-genes_stages_I[genes_stages_I %in% genes_Stage_I$gene]
genes_stages_II<-genes_stages_II[genes_stages_II %in% genes_Stage_II$gene]
genes_stages_III<-genes_stages_III[genes_stages_III %in% genes_Stage_III$gene]

p_stage_I_unpaired<-ggplot(unstranded_data_samples_unapaired[unstranded_data_samples_unapaired$gene_id %in% genes_stages_I,], aes(x=stages, y=RPKM, fill=stages)) +   geom_boxplot()+ facet_wrap(~SYMBOL.x, nrow = 5,ncol = 6, scales="free")+ theme_bw() + ggtitle("Biomarkers for stage I") + stat_compare_means(comparisons = my_comparisons, vjust =0.03, label.x.npc="bottom", method="t.test") + scale_fill_manual(values=c('#e8f2a1', '#729fcf', '#ffaa95')) + stat_summary(fun.y=mean, geom="point", shape=20, size=8, color="red", fill="red")+ geom_jitter(color="black", size=0.4, alpha=0.9)
p_stage_II_unpaired<-ggplot(unstranded_data_samples_unapaired[unstranded_data_samples_unapaired$gene_id %in% genes_stages_II,], aes(x=stages, y=RPKM, fill=stages)) +   geom_boxplot()+ facet_wrap(~SYMBOL.x, nrow = 3,ncol = 3, scales="free")+ theme_bw() + ggtitle("Biomarkers for stage II") + stat_compare_means(comparisons = my_comparisons, vjust =0.03, label.x.npc="bottom", method="t.test") + scale_fill_manual(values=c('#e8f2a1', '#729fcf', '#ffaa95'))+ stat_summary(fun.y=mean, geom="point", shape=20, size=8, color="red", fill="red")+ geom_jitter(color="black", size=0.4, alpha=0.9)
p_stage_III_unpaired<-ggplot(unstranded_data_samples_unapaired[unstranded_data_samples_unapaired$gene_id %in% genes_stages_III,], aes(x=stages, y=RPKM, fill=stages)) +   geom_boxplot()+ facet_wrap(~SYMBOL.x, nrow = 3,ncol = 3, scales="free")+ theme_bw() + ggtitle("Biomarkers for stage III") + stat_compare_means(comparisons = my_comparisons, vjust =0.03, label.x.npc="bottom", method="t.test") + scale_fill_manual(values=c('#e8f2a1', '#729fcf', '#ffaa95'))+ stat_summary(fun.y=mean, geom="point", shape=20, size=8, color="red", fill="red") + geom_jitter(color="black", size=0.4, alpha=0.9)
