library(pathview)

#map04012 ErbB signaling pathway			04012
#N01579 CD80/CD86-CTLA4-PP2A signaling pathway		04660
#map04668 TNF signaling pathway				04668
#map04151 PI3K-Akt signaling pathway			04151
#map04115 p53 signaling pathway				04115
#map04010 MAPK signaling pathway			04010
#map04630 JAK-STAT signaling pathway			04630
#nt06507 TGFB signaling					06507
#N00151 TNF-NFKB signaling pathway			04668
#map04330 Notch signaling pathway			04330
#map04340 Hedgehog signaling pathway			04340
#map01521 EGFR tyrosine kinase inhibitor resistance	01521	
#map04350 TGF-beta signaling pathway			04350
#map04370 VEGF signaling pathway			04370
#map03320 PPAR signaling pathway			03320
#map04310 Wnt signaling pathway				04310

mypathway_ErbB<-"04012"
mypathway_CTLA4<-"04660"
mypathway_TNF<-"04668"
mypathway_PI3K<-"04151"
mypathway_MAPK<-"04010"
mypathway_JAK<-"04630"
mypathway_TGFB<-"06507"
mypathway_TNF<-"04668"
mypathway_Notch<-"04330"
mypathway_Hedgehog<-"04340"
mypathway_EGFR<-"01521"
mypathway_TGF<- "04350"
mypathway_VEGF<-"04370"
mypathway_WNT<-"04310"
mypathway_PPAR<-"04320"



genes<-c(genes_unique_Stage_I$gene_id,genes_unique_Stage_II$gene_id,genes_unique_Stage_III$gene_id)

# For each gene, add gene_id
for (gene_row in rownames(log2change_tumor_control))
{	
	# Store gene id in the vector
	# Simply trim the gene id before the "." to save it in the ENSEML format
	log2change_tumor_control[gene_row,"ENSEMBL"]<-strsplit(log2change_tumor_control[gene_row,"gene"], split = "\\.")[[1]][1]	
}

logFC<-log2change_tumor_control[log2change_tumor_control$ENSEMBL %in% as.vector(genes),"log2change"]
ids_stage_I      <-bitr(log2change_tumor_control$ENSEMBL, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
names(logFC)<-ids_stage_I[ids_stage_I$ENSEMBL %in% genes,"ENTREZID"]


mypathway_ErbB<-"04012"
mypathway_CTLA4<-"04660"
mypathway_TNF<-"04668"
mypathway_PI3K<-"04151"
mypathway_MAPK<-"04010"
mypathway_JAK<-"04630"
mypathway_TGFB<-"06507"
mypathway_TNF<-"04668"
mypathway_Notch<-"04330"
mypathway_Hedgehog<-"04340"
mypathway_EGFR<-"01521"
mypathway_TGF<- "04350"
mypathway_VEGF<-"04370"

pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_ErbB)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_CTLA4)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_TNF)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_PI3K)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_MAPK)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_JAK)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_TGFB)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_TNF)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_Notch)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_Hedgehog)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_EGFR)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_TGF)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_VEGF)
# Reactome WNT signalling pathway is enriched with the stage-specific genes for stage III.
# Choosen because it is abundant in numeber of stage-specfic genes (13 genes): DVL2, DVL2, PSMB1, AP2S1, PSMC6, PSMA7, PSMD7, AKT2, BCL9, CLTA, PSMD4, PSMD4, PSMC2
# Kegg pathway gene shows genes DVL1, FRP, BAMBI biomarkers acting in this patwhay.
# FRP appears repressing the activation of Wnt to Frizzled-LRP5/6
# ids_stage_I[ids_stage_I$ENSEMBL =="ENSG00000107404",] # DVL1
# ids_stage_I[ids_stage_I$ENSEMBL =="ENSG00000104332",] # FRP
# ids_stage_I[ids_stage_I$ENSEMBL =="ENSG00000095739",] # BAMBI
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_WNT)
pathview(gene.data=logFC*100,species="hsa",pathway="03320")


#######################################################################################################################################
# Path to files of selected_genes                                                                                                             # 
selected_genes_Stage_I_file       <-paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_Stage_","sample_stage_I",".tsv",sep="")
selected_genes_Stage_II_file      <-paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_Stage_","sample_stage_II",".tsv",sep="")
selected_genes_Stage_III_file     <-paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_Stage_","sample_stage_III",".tsv",sep="")
#######################################################################################################################################
# Load data                                                                                                                           #
selected_genes_Stage_I_data       <-read.table(file = selected_genes_Stage_I_file, sep = '\t', header = TRUE,fill=TRUE)                        #
selected_genes_Stage_II_data      <-read.table(file = selected_genes_Stage_II_file, sep = '\t', header = TRUE,fill=TRUE)                       #
selected_genes_Stage_III_data     <-read.table(file = selected_genes_Stage_III_file, sep = '\t', header = TRUE,fill=TRUE)                      #
                                                                                                                                      #
# Set rownames                                                                                                                        #
rownames(selected_genes_Stage_I_data)<-selected_genes_Stage_I_data$gene                                                               #
rownames(selected_genes_Stage_II_data)<-selected_genes_Stage_II_data$gene                                                             #
rownames(selected_genes_Stage_III_data)<-selected_genes_Stage_III_data$gene                                                           #
#######################################################################################################################################
selected_genes_Stage_I_data$Stage<-"Stage I"
selected_genes_Stage_II_data$Stage<-"Stage II"
selected_genes_Stage_III_data$Stage<-"Stage III"
selected_genes_Stage_data<-rbind(selected_genes_Stage_I_data,selected_genes_Stage_II_data, selected_genes_Stage_III_data)
selected_genes_Stage_data$ENSEMBL<-""

# For each gene, add gene_id
for (gene_row in rownames(selected_genes_Stage_data))
{	
	# Store gene id in the vector
	# Simply trim the gene id before the "." to save it in the ENSEML format
	selected_genes_Stage_data[gene_row,"ENSEMBL"]<-strsplit(selected_genes_Stage_data[gene_row,"gene"], split = "\\.")[[1]][1]
}

# For each gene, add gene_id
for (gene_row in genes)
{	
	# Store gene id in the vector
	# Simply trim the gene id before the "." to save it in the ENSEML format
	genes_ids_all<-c(genes_ids_all,strsplit(gene_row, split = "\\.")[[1]][1])
}

genes<-unique(c(rownames(selected_genes_Stage_I_data),rownames(selected_genes_Stage_II_data),rownames(selected_genes_Stage_III_data)))
genes_ids_all<-c()

# For each gene, add gene_id
for (gene_row in genes)
{	
	# Store gene id in the vector
	# Simply trim the gene id before the "." to save it in the ENSEML format
	genes_ids_all<-c(genes_ids_all,strsplit(gene_row, split = "\\.")[[1]][1])
}
log2change_tumor_control$ENSEMBL<-""
# For each gene, add gene_id
for (gene_row in rownames(log2change_tumor_control))
{	
	# Store gene id in the vector
	# Simply trim the gene id before the "." to save it in the ENSEML format
	log2change_tumor_control[gene_row,"ENSEMBL"]<-strsplit(log2change_tumor_control[gene_row,"gene"], split = "\\.")[[1]][1]
}
logFC<-log2change_tumor_control[log2change_tumor_control$ENSEMBL %in% genes_ids_all,"log2change"]
ids_stage_I      <-bitr(log2change_tumor_control$ENSEMBL, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
names(logFC)<-ids_stage_I[ids_stage_I$ENSEMBL %in% genes_ids_all,"ENTREZID"]


pathview(gene.data=logFC*100,species="hsa",pathway="05223")
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_ErbB)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_CTLA4)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_TNF)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_PI3K)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_MAPK)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_JAK)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_TGFB)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_TNF)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_Notch)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_Hedgehog)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_EGFR)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_TGF)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_VEGF)
pathview(gene.data=logFC*100,species="hsa",pathway="03320")


logFC_1<-logFC[ids_stage_I[ids_stage_I$SYMBOL %in% c("DVL2","PHB1","ELOC", "PSMC6","SEH1L"),"ENTREZID"]]*100
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_WNT)

