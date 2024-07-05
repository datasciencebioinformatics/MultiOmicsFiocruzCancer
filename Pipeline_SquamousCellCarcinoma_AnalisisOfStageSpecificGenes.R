# Reat stage specific genes 
stage_specific_genes <-  read.xlsx(file="/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Table3.xlsx", 2)   # read first sheet

# Genes that are tumor genes (RPKM ≥ 4, tumor vs. normal samples : log2foldchange >1 and paired t-test FDR ≤ 0.05) but also stage-specific genes (RPKM ≥ 4, stage-specific vs. normal samples : log2foldchange >1 and paired t-test FDR ≤ 0.05).

# (RPKM ≥ 4, tumor vs. normal samples : log2foldchange >1 and paired t-test FDR ≤ 0.05)
tumor_genes<-stage_specific_genes[stage_specific_genes$tumor.normal.log2fc >1 & stage_specific_genes$tumor.normal.FDR < 0.05,]

# Stage specific genes
stage_specific_genes_stage_I   <-stage_specific_genes[stage_specific_genes$Stage_I.normal.log2fc >1 & stage_specific_genes$Stage_I.normal.FDR < 0.05,]
stage_specific_genes_stage_II  <-stage_specific_genes[stage_specific_genes$Stage_II.normal.log2fc >1 & stage_specific_genes$Stage_II.normal.FDR < 0.05,]
stage_specific_genes_stage_III <-stage_specific_genes[stage_specific_genes$Stage_III.normal.log2fc >1 & stage_specific_genes$Stage_III.normal.FDR < 0.05,]
 

# Stage I, comparisson statistics against Stage II and Stage III
# Stage I   :  Gene is selected if Stage I > Stage II and  Stage I > Stage III
# Stage II  :  Gene is selected if Stage II > Stage I and  Stage II > Stage III
# Stage III : selected if Stage III > Stage I and  Stage III > Stage II

# Stage specific genes
stage_I_selection<-head(stage_specific_genes[order(-stage_specific_genes$Stage_I.normal.log2fc),],n=20)
stage_II_selection<-head(stage_specific_genes[order(-stage_specific_genes$Stage_II.normal.log2fc),],n=20)
stage_III_selection<-head(stage_specific_genes[order(-stage_specific_genes$Stage_III.normal.log2fc),],n=20)


p_stage_I_unpaired<-ggplot(unstranded_data_samples_unapaired[unstranded_data_samples_unapaired$ENSEMBL.x %in% stage_I_selection$ENSEMBL,], aes(x=stages, y=RPKM, fill=stages)) +   geom_boxplot()+ facet_wrap(~SYMBOL.x, nrow = 5,ncol = 6, scales="free")+ theme_bw() + ggtitle("Biomarkers for stage I") + stat_compare_means(comparisons = my_comparisons, vjust =0.03, label.x.npc="bottom", method="t.test") + scale_fill_manual(values=c('#e8f2a1', '#729fcf', '#ffaa95'))
p_stage_II_unpaired<-ggplot(unstranded_data_samples_unapaired[unstranded_data_samples_unapaired$ENSEMBL.x %in% stage_II_selection$ENSEMBL,], aes(x=stages, y=RPKM, fill=stages)) +   geom_boxplot()+ facet_wrap(~SYMBOL.x, nrow = 3,ncol = 3, scales="free")+ theme_bw() + ggtitle("Biomarkers for stage II") + stat_compare_means(comparisons = my_comparisons, vjust =0.03, label.x.npc="bottom", method="t.test") + scale_fill_manual(values=c('#e8f2a1', '#729fcf', '#ffaa95'))
p_stage_III_unpaired<-ggplot(unstranded_data_samples_unapaired[unstranded_data_samples_unapaired$ENSEMBL.x %in% stage_III_selection$ENSEMBL,], aes(x=stages, y=RPKM, fill=stages)) +   geom_boxplot()+ facet_wrap(~SYMBOL.x, nrow = 3,ncol = 3, scales="free")+ theme_bw() + ggtitle("Biomarkers for stage II") + stat_compare_means(comparisons = my_comparisons, vjust =0.03, label.x.npc="bottom", method="t.test") + scale_fill_manual(values=c('#e8f2a1', '#729fcf', '#ffaa95'))
