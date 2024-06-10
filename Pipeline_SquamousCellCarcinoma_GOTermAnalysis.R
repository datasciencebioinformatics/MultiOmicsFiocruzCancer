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
go_ALL_Stage_I    <- enrichGO(gene = ids_stage_I$SYMBOL,  OrgDb  = org.Hs.eg.db,      ont = "ALL", pAdjustMethod = "BH",pvalueCutoff  = 0.1,qvalueCutoff  = 0.15,readable = TRUE, minGSSize = 3,keyType = "SYMBOL")
go_ALL_Stage_II   <- enrichGO(gene = ids_stage_II$SYMBOL, OrgDb  = org.Hs.eg.db,      ont = "ALL",  pAdjustMethod = "BH",pvalueCutoff  = 0.1,qvalueCutoff  = 0.15,readable = TRUE, minGSSize = 3,keyType = "SYMBOL")
go_ALL_Stage_III  <- enrichGO(gene = ids_stage_III$SYMBOL,OrgDb  = org.Hs.eg.db,      ont = "ALL",   pAdjustMethod = "BH",pvalueCutoff  = 0.1,qvalueCutoff  = 0.15,readable = TRUE, minGSSize = 3,keyType = "SYMBOL")
go_ALL_Stage      <- enrichGO(gene = genes_Stage_ALL$SYMBOL,OrgDb  = org.Hs.eg.db,    ont = "ALL",   pAdjustMethod = "BH",pvalueCutoff  = 0.1,qvalueCutoff  = 0.15,readable = TRUE, minGSSize = 3,keyType = "SYMBOL")
########################################################################################################################################
kegg_ALL_Stage_I    <- enrichKEGG(gene = ids_stage_I$ENTREZ,  organism     = 'hsa',    pvalueCutoff = 0.15)
kegg_ALL_Stage_II   <- enrichKEGG(gene = ids_stage_II$ENTREZ,  organism     = 'hsa',    pvalueCutoff = 0.15)
kegg_ALL_Stage_III  <- enrichKEGG(gene = ids_stage_III$ENTREZ,  organism     = 'hsa',    pvalueCutoff = 0.15)
kegg_ALL_Stage      <- enrichKEGG(gene = genes_Stage_ALL$ENTREZ,  organism     = 'hsa',    pvalueCutoff = 0.15)
########################################################################################################################################
reactome_ALL_Stage_I    <- enrichPathway(gene=ids_stage_I$ENTREZ, pvalueCutoff = 0.15, readable=TRUE)
reactome_ALL_Stage_II   <- enrichPathway(gene=ids_stage_II$ENTREZ, pvalueCutoff = 0.15, readable=TRUE)
reactome_ALL_Stage_III  <- enrichPathway(gene=ids_stage_III$ENTREZ, pvalueCutoff = 0.15, readable=TRUE)
reactome_ALL_Stage      <- enrichPathway(gene = genes_Stage_ALL$ENTREZ,pvalueCutoff = 0.05, readable=TRUE)
########################################################################################################################################
# A table with the associate go term per gene.
# First, take the results per stage
# Go terms
go_results_Stage_I<-go_ALL_Stage_I@result
go_results_Stage_II<-go_ALL_Stage_II@result
go_results_Stage_III<-go_ALL_Stage_III@result

# Kegg terms
kegg_results_Stage_I<-kegg_ALL_Stage_I@result
kegg_results_Stage_II<-kegg_ALL_Stage_II@result
kegg_results_Stage_III<-kegg_ALL_Stage_III@result

# Reactome terms
reactome_results_Stage_I<-reactome_ALL_Stage_I@result
reactome_results_Stage_II<-reactome_ALL_Stage_II@result
reactome_results_Stage_III<-reactome_ALL_Stage_III@result

# Second, add information about stage
go_results_Stage_I$Stage<-"Stage I"
go_results_Stage_II$Stage<-"Stage II"
go_results_Stage_III$Stage<-"Stage III"

# Second, add information about stage
kegg_results_Stage_I$Stage<-"Stage I"
kegg_results_Stage_II$Stage<-"Stage II"
kegg_results_Stage_III$Stage<-"Stage III"

# Second, add information about stage
reactome_results_Stage_I$Stage<-"Stage I"
reactome_results_Stage_II$Stage<-"Stage II"
reactome_results_Stage_III$Stage<-"Stage III"

# Merge all stages
go_results_all_Stages         <-rbind(go_results_Stage_I,go_results_Stage_II,go_results_Stage_III)[,c("ID","p.adjust","Description","geneID","Count","Stage")]
kegg_results_all_Stages       <-rbind(kegg_results_Stage_I,kegg_results_Stage_II,kegg_results_Stage_III)[,c("ID","p.adjust","Description","geneID","Count","Stage")]
reactome_results_all_Stages   <-rbind(reactome_results_Stage_I,reactome_results_Stage_II,reactome_results_Stage_III)[,c("ID","p.adjust","Description","geneID","Count","Stage")]

# Merge all stages
go_results_all_Stages$Layer<-"GO"
kegg_results_all_Stages$Layer<-"KEGG"
reactome_results_all_Stages$Layer<-"Reactome"

# Merge all tables
all_anotation_results<-rbind(go_results_all_Stages,kegg_results_all_Stages,reactome_results_all_Stages)

# Filter by padj
#all_anotation_results<-all_anotation_results[which(all_anotation_results$p.adjust<0.15),]
######################################################################################################################
# A data.frame to store all results
df_all_annotation<-data.frame(ID=c(),p.adjust=c(),Description=c(),geneID=c(), SYMBOL=c(),Count=c(),Stage=c(),Layer=c())

# For each annotation
for (annotation in rownames(all_anotation_results))
{

	# Save all variables
	ID            <- all_anotation_results[annotation,"ID"]   
	p.adjust      <- all_anotation_results[annotation,"p.adjust"]   
	Description   <- all_anotation_results[annotation,"Description"]   
	geneID        <- all_anotation_results[annotation,"geneID"]   
	Count         <- all_anotation_results[annotation,"Count"]   
	Stage         <- all_anotation_results[annotation,"Stage"]   
	Layer         <- all_anotation_results[annotation,"Layer"]   

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
# Save file 
write.xlsx(x=df_all_annotation,file=paste(output_dir,"df_all_annotation",".xlsx",sep=""), sheet="annnotation")
######################################################################################################################
# Load interactome
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadInteractomeUniqueGenes.R")

# Load interactome
genes_in_interactome<-unique(genes_Stage_ALL[genes_Stage_ALL$gene_id %in% intersect(genes_Stage_ALL$gene_id,c(interactome_all_stage$Gene1,interactome_all_stage$Gene2)),])

# Set rownames as gene_ids
rownames(genes_in_interactome)<-genes_in_interactome$gene_id

# Set gene symbols
interactome_all_stage$Gene1<-genes_in_interactome[interactome_all_stage$Gene1,"SYMBOL"]
interactome_all_stage$Gene2<-genes_in_interactome[interactome_all_stage$Gene2,"SYMBOL"]

# gene annotation
interactome_annotation<-data.frame(ID=rownames(interactome_all_stage),p.adjust=0,Description=interactome_all_stage$Gene1,geneID=interactome_all_stage$Gene2, SYMBOL=interactome_all_stage$Gene2,Count=0,Stage=interactome_all_stage$Stage,Layer="interactome")
######################################################################################################################
# Stage II
Stage="Stage III"

# Store annotation of stage II
df_all_annotation_per_stage<-df_all_annotation[which(df_all_annotation$Stage==Stage),]

# interactome_annotation
interactome_annotation_stage<-interactome_annotation[which(interactome_annotation$Stage==Stage),]

# interactome_annotation_stage
selected_interactome<-c(paste(interactome_annotation_stage$Description,interactome_annotation_stage$SYMBOL,sep="-"),
paste(interactome_annotation_stage$SYMBOL,interactome_annotation_stage$Description,sep="-"))
######################################################################################################################
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

selection_GO          <-names(tail(sort(table_GO),n=10))
selection_KEGG        <-names(tail(sort(table_KEGG),n=10))
selection_REACTOM     <-names(tail(sort(table_REACTOME),n=10))

df_all_annotation_selected_pathways<-rbind(df_all_annotation[which(df_all_annotation$Description %in% selection_GO),],
df_all_annotation[which(df_all_annotation$Description %in% selection_KEGG),],
df_all_annotation[which(df_all_annotation$Description %in% selection_REACTOM),],
interactome_annotation_stage)
####################################################################################################################
# All stages
graph_all_stages <- graph_from_data_frame(d=df_all_annotation_selected_pathways[,c("SYMBOL","Description")], vertices=unique(c(df_all_annotation_selected_pathways$SYMBOL,df_all_annotation_selected_pathways$Description)), directed=F)  

# df_all_eges
df_all_eges<-data.frame(get.edgelist(graph_all_stages))

# Set ronames
df_all_eges$names<-paste(df_all_eges$X1,df_all_eges$X2,sep="-")
####################################################################################################################
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

# Set size of the node according to the dregree
V(graph_all_stages)$size                                                                           <- degree(graph_all_stages)

# Vertice colours of genes
E(graph_all_stages)$color                                                                         <- "lightgray"

# Set ronames
E(graph_all_stages)[which(df_all_eges$names %in% selected_interactome)]$color                     <- "black"

# FindClusters_resolution
png(filename=paste(output_folder,"Plot_Stage_III.png",sep=""), width = 30, height = 30, res=600, units = "cm")
	#plot(graph_all_stages, layout=layout_nicely) # Stage I
	plot(graph_all_stages, layout=  layout_with_kk) # Stage II
	#plot(graph_all_stages, layout=   layout_nicely) # Stage II
dev.off()

# Save file 
write.xlsx(x=df_all_annotation_selected_pathways,file=paste(output_dir,"df_all_annotation",".xlsx",sep=""), sheet="Stage III")


######################################################################################################################
GO=
KEGG=
REACTOME=
Gene=

png(filename=paste(output_folder,"legend.png",sep=""), width = 5, height = 5, res=600, units = "cm")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c('GO', 'KEGG', 'REACTOME','Gene'), pch=16, pt.cex=3, cex=1.5, bty='n',col = c('#ff6347', '#ffd700', '#7f7f7f', '#0072b2'))
mtext("Legend", at=0.2, cex=2)
dev.off()
