############################################################################################################################################
# cyrestGET(operation = NULL, parameters = NULL, base.url = "http://localhost:1234")
# Analyses with the combination of parameter line 119 of the document Parametrization.xlsx
# ≥3	≥1	≤0.05	≥0.85	5456	1798/25	1887/70	1991/182	204/191/1.3396	225/207/1.4054	242/206/1.2978	1276/3819/3.7205	1345/4143/3.7816	1440/4646/3.8299
# ENSEMBL ids were converted to ENTREZ ids. enrichGO on org.Hs.eg.db was used (pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05, minGSSize = 3) to anotate 23, 62 and 169 genes from stages I, II and III, respectivelly. Then cnetplot was used to show asociations of genes to top 10 categories.
#######################################################################################################################################
file_unique_gene_stages_I    <- paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_I",".tsv",sep="")
file_unique_gene_stages_II   <- paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_II",".tsv",sep="")
file_unique_gene_stages_III  <- paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_III",".tsv",sep="")
#######################################################################################################################################
genes_unique_Stage_I         <-read.table(file = file_unique_gene_stages_I, sep = '\t', header = TRUE,fill=TRUE)    #
genes_unique_Stage_II        <-read.table(file = file_unique_gene_stages_II, sep = '\t', header = TRUE,fill=TRUE)    #
genes_unique_Stage_III       <-read.table(file = file_unique_gene_stages_III, sep = '\t', header = TRUE,fill=TRUE)    #
#######################################################################################################################################
# gene_id
genes_Stage_I$gene_id                   <-""
genes_Stage_II$gene_id                  <-""
genes_Stage_III$gene_id                 <-""
#######################################################################################################################################
# Calculate gse all stages
# For each gene in stage I
genes_ids_stage_I<-c()
genes_ids_stage_II<-c()
genes_ids_stage_III<-c()
genes_ids_all<-c()

# For each gene, add gene_id
for (gene_row in rownames(genes_unique_Stage_I))
{	
	# Store gene id in the vector
	genes_unique_Stage_I[gene_row,"gene_id"]<-strsplit(genes_unique_Stage_I[gene_row,"gene"], split = "\\.")[[1]][1]
	
}
# For each gene, add gene_id
for (gene_row in rownames(genes_unique_Stage_II))
{	
	# Store gene id in the vector
	genes_unique_Stage_II[gene_row,"gene_id"]<-strsplit(genes_unique_Stage_II[gene_row,"gene"], split = "\\.")[[1]][1]	
}
# For each gene, add gene_id
for (gene_row in rownames(genes_unique_Stage_III))
{	
	# Store gene id in the vector
	genes_unique_Stage_III[gene_row,"gene_id"]<-strsplit(genes_unique_Stage_III[gene_row,"gene"], split = "\\.")[[1]][1]	
}
########################################################################################################################################
# ids_stage_I
ids_stage_I    <-bitr(genes_unique_Stage_I$gene_id, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
ids_stage_II   <-bitr(genes_unique_Stage_II$gene_id, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
ids_stage_III  <-bitr(genes_unique_Stage_III$gene_id, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
########################################################################################################################################
colnames(ids_stage_I)   <-c("gene_id","ENTREZID","SYMBOL")
colnames(ids_stage_II)  <-c("gene_id","ENTREZID","SYMBOL")
colnames(ids_stage_III) <-c("gene_id","ENTREZID","SYMBOL")
########################################################################################################################################
genes_Stage_I  <-merge(genes_unique_Stage_I,ids_stage_I,by="gene_id")
genes_Stage_II <-merge(genes_unique_Stage_II,ids_stage_II,by="gene_id")
genes_Stage_III<-merge(genes_unique_Stage_III,ids_stage_III,by="gene_id")
genes_Stage_ALL<-unique(rbind(ids_stage_I,ids_stage_II,ids_stage_III))
rownames(genes_Stage_ALL)<-genes_Stage_ALL$ENTREZID
########################################################################################################################################
go_ALL_Stage            <- enrichGO(gene = genes_Stage_ALL$SYMBOL,OrgDb  = org.Hs.eg.db,    ont = "ALL",   pAdjustMethod = "BH",pvalueCutoff  = 0.1,qvalueCutoff  = 0.15,readable = TRUE, minGSSize = 3,keyType = "SYMBOL")
kegg_ALL_Stage          <- enrichKEGG(gene = genes_Stage_ALL$ENTREZ,  organism     = 'hsa',    pvalueCutoff = 0.15)
reactome_ALL_Stage      <- enrichPathway(gene = genes_Stage_ALL$ENTREZ,pvalueCutoff = 0.05, readable=TRUE)
########################################################################################################################################
# A table with the associate go term per gene.
# GO, KEGG e REACTOME
go_ALL_Stage       <-go_ALL_Stage@result[,c("ID","p.adjust","Description","geneID","Count")]
kegg_ALL_Stage     <-kegg_ALL_Stage@result[,c("ID","p.adjust","Description","geneID","Count")]
reactome_ALL_Stage <-reactome_ALL_Stage@result[,c("ID","p.adjust","Description","geneID","Count")]

# Merge all stages
go_ALL_Stage$Layer       <-"GO"
kegg_ALL_Stage$Layer     <-"KEGG"
reactome_ALL_Stage$Layer <-"Reactome"

# Merge all stages
annotation_all_layers<-rbind(go_ALL_Stage,kegg_ALL_Stage,reactome_ALL_Stage)

# Filter by padj
annotation_all_layers<-annotation_all_layers[which(annotation_all_layers$p.adjust<0.25),]
######################################################################################################################
# A data.frame to store all results
df_all_annotation<-data.frame(ID=c(),p.adjust=c(),Description=c(),geneID=c(), SYMBOL=c(),Count=c(),Stage=c(),Layer=c())

# For each annotation
for (annotation_row in  1:dim(annotation_all_layers)[1])
{
	# Save all variables
	ID            <- all_anotation_results[annotation_row,"ID"]   
	p.adjust      <- all_anotation_results[annotation_row,"p.adjust"]   
	Description   <- all_anotation_results[annotation_row,"Description"]   
	geneID        <- all_anotation_results[annotation_row,"geneID"]   
	Count         <- all_anotation_results[annotation_row,"Count"]   
	Stage         <- "ALL"
	Layer         <- all_anotation_results[annotation_row,"Layer"]   

	# All genes of stage
	all_genes_of_annotation<-strsplit(x=geneID,"/",fixed=T)[[1]]

	# For all genes
	for (genes in all_genes_of_annotation)
	{
		# If layer equal to Kegg
		if(Layer=="KEGG")
		{
			# Make the id converstion			
			gene_symbol<-genes_Stage_ALL[which(genes_Stage_ALL$ENTREZID %in% genes),"SYMBOL"]
			genes_id<-genes_Stage_ALL[which(genes_Stage_ALL$ENTREZID %in% genes),"gene_id"]
		}else # If not kegg 
		{
			# Make the id converstion			
			gene_symbol<-genes_Stage_ALL[which(genes_Stage_ALL$SYMBOL %in% genes),"SYMBOL"]
			genes_id<-genes_Stage_ALL[which(genes_Stage_ALL$SYMBOL %in% genes),"gene_id"]
		}
		# gene annotation
		gene_annotation<-data.frame(ID=ID,p.adjust=p.adjust,Description=Description,geneID=genes_id, SYMBOL=gene_symbol,Count=Count,Stage=Stage,Layer=Layer)

		# df_all_annotation
		df_all_annotation<-rbind(df_all_annotation,gene_annotation)
	}	
}
######################################################################################################################
# Set all stages
df_all_annotation$Stage<-""
df_all_annotation[which(df_all_annotation$geneID %in% genes_unique_Stage_I$gene_id),"Stage"]<-"Stage I"
df_all_annotation[which(df_all_annotation$geneID %in% genes_unique_Stage_II$gene_id),"Stage"]<-"Stage II"
df_all_annotation[which(df_all_annotation$geneID %in% genes_unique_Stage_III$gene_id),"Stage"]<-"Stage III"

# Save file 
write.xlsx(x=df_all_annotation,file=paste(output_dir,"df_all_annotation_V2",".xlsx",sep=""), sheet="annnotation")
######################################################################################################################
# Load interactome
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadInteractomeUniqueGenes.R")
######################################################################################################################
# Stage I
df_all_annotation_per_stage<-df_all_annotation[which(df_all_annotation$Stage=="Stage I"),]

# A table for the count of go terms
table_GO<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="GO","Description"])
table_KEGG<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="KEGG","Description"])
table_REACTOME<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="Reactome","Description"])

# A vector with all gene symbols
gene_symbols<-df_all_annotation_per_stage$SYMBOL

# A vector with all gene symbols
annotation_description_GO       <-unique(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="GO","Description"])
annotation_description_KEGG     <-unique(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="KEGG","Description"])
annotation_description_REACTOME <-unique(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="Reactome","Description"])

# A vector with all gene symbols
gene_Stage_I                 <-unique(df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage I"),"SYMBOL"])
gene_Stage_II                <-unique(df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage II"),"SYMBOL"])
gene_Stage_III               <-unique(df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage III"),"SYMBOL"])

# All stages
graph_all_stages <- graph_from_data_frame(d=df_all_annotation_per_stage[,c("SYMBOL","Description")], vertices=unique(c(df_all_annotation_per_stage$SYMBOL,df_all_annotation_per_stage$Description)), directed=F) 


# Vertice colours of genes
V(graph_all_stages)$color                                                                         <- "#0072b2"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_GO)]$color       <- "#ff6347"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_KEGG)]$color     <- "#ffd700"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_REACTOME)]$color <- "#7f7f7f"

# Vertice colours of genes
V(graph_all_stages)$shape                                                                           <-"circle"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_GO)]$shape        <- "square"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_KEGG)]$shape      <- "square"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_REACTOME)]$shape  <- "square"

# FindClusters_resolution
png(filename=paste(output_folder,"Plot_Stage_I.png",sep=""), width = 25, height = 25, res=600, units = "cm")
	plot(graph_all_stages, layout=layout_with_kk)
dev.off()










######################################################################################################################
# Stage I
df_all_annotation_per_stage<-df_all_annotation[which(df_all_annotation$Stage=="Stage II"),]



# A table for the count of go terms
table_GO<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="GO","Description"])
table_KEGG<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="KEGG","Description"])
table_REACTOME<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="Reactome","Description"])

# A vector with all gene symbols
gene_symbols<-df_all_annotation_per_stage$SYMBOL

# A vector with all gene symbols
annotation_description_GO       <-unique(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="GO","Description"])
annotation_description_KEGG     <-unique(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="KEGG","Description"])
annotation_description_REACTOME <-unique(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="Reactome","Description"])

# A vector with all gene symbols
gene_Stage_I                 <-unique(df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage I"),"SYMBOL"])
gene_Stage_II                <-unique(df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage II"),"SYMBOL"])
gene_Stage_III               <-unique(df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage III"),"SYMBOL"])

# All stages
graph_all_stages <- graph_from_data_frame(d=df_all_annotation_per_stage[,c("SYMBOL","Description")], vertices=unique(c(df_all_annotation_per_stage$SYMBOL,df_all_annotation_per_stage$Description)), directed=F) 


# Vertice colours of genes
V(graph_all_stages)$color                                                                         <- "#0072b2"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_GO)]$color       <- "#ff6347"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_KEGG)]$color     <- "#ffd700"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_REACTOME)]$color <- "#7f7f7f"

# Vertice colours of genes
V(graph_all_stages)$shape                                                                           <-"circle"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_GO)]$shape        <- "square"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_KEGG)]$shape      <- "square"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_REACTOME)]$shape  <- "square"

# FindClusters_resolution
png(filename=paste(output_folder,"Plot_Stage_II.png",sep=""), width = 25, height = 25, res=1200, units = "cm")
	plot(graph_all_stages, layout=layout_with_kk)
dev.off()




######################################################################################################################
# Stage I
df_all_annotation_per_stage<-df_all_annotation[which(df_all_annotation$Stage=="Stage III"),]

# A table for the count of go terms
table_GO<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="GO","Description"])
table_KEGG<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="KEGG","Description"])
table_REACTOME<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="Reactome","Description"])

# A vector with all gene symbols
gene_symbols<-df_all_annotation_per_stage$SYMBOL

# A vector with all gene symbols
annotation_description_GO       <-unique(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="GO","Description"])
annotation_description_KEGG     <-unique(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="KEGG","Description"])
annotation_description_REACTOME <-unique(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="Reactome","Description"])

# A vector with all gene symbols
gene_Stage_I                 <-unique(df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage I"),"SYMBOL"])
gene_Stage_II                <-unique(df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage II"),"SYMBOL"])
gene_Stage_III               <-unique(df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage III"),"SYMBOL"])

# All stages
graph_all_stages <- graph_from_data_frame(d=df_all_annotation_per_stage[,c("SYMBOL","Description")], vertices=unique(c(df_all_annotation_per_stage$SYMBOL,df_all_annotation_per_stage$Description)), directed=F) 


# Vertice colours of genes
V(graph_all_stages)$color                                                                         <- "#0072b2"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_GO)]$color       <- "#ff6347"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_KEGG)]$color     <- "#ffd700"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_REACTOME)]$color <- "#7f7f7f"

# Vertice colours of genes
V(graph_all_stages)$shape                                                                           <-"circle"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_GO)]$shape        <- "square"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_KEGG)]$shape      <- "square"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_REACTOME)]$shape  <- "square"

# FindClusters_resolution
png(filename=paste(output_folder,"Plot_Stage_III.png",sep=""), width = 25, height = 25, res=600, units = "cm")
	plot(graph_all_stages, layout=layout_with_kk)
dev.off()
