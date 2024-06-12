annotation_stages_all
######################################################################################################################
# A data.frame to store all results
df_all_annotation<-data.frame(gene_id=c(),gene=c(),log2change=c(),Category=c(),Pvalue=c(),FDR=c(),ENTREZID=c(),Symbol=c(),Description=c(),genecards_Category=c(),UniProt_ID=c(),GIFtS=c(),GC_id=c(),GeneCards_Summary=c(),Stage=c(),CluterProfiler=c() )

# For each annotation
for (annotation in rownames(annotation_stages_all))
{

	# Save all variables
	gene_id      <- annotation_stages_all[annotation,"gene_id"]   
	gene         <- annotation_stages_all[annotation,"gene"]   
	log2change   <- annotation_stages_all[annotation,"log2change"]   
	Category     <- annotation_stages_all[annotation,"Category"]   
	Pvalue       <- annotation_stages_all[annotation,"Pvalue"]   
	FDR          <- annotation_stages_all[annotation,"FDR"]   
	ENTREZID     <- annotation_stages_all[annotation,"ENTREZID"]   
	Symbol       <- annotation_stages_all[annotation,"Symbol"]   
	Description  <- annotation_stages_all[annotation,"Description"]   
	genecards_Category   <- annotation_stages_all[annotation,"genecards_Category"]   
	UniProt_ID           <- annotation_stages_all[annotation,"UniProt_ID"]   
	GIFtS                <- annotation_stages_all[annotation,"GIFtS"]   
	GC_id                <- annotation_stages_all[annotation,"GC_id"]   
	GeneCards_Summary    <- annotation_stages_all[annotation,"GeneCards_Summary"]   
	Stage                <- annotation_stages_all[annotation,"Stage"]   
	CluterProfiler	     <- annotation_stages_all[annotation,"CluterProfiler"]   

	# All genes of stage
	all_genes_of_annotation<-strsplit(x=CluterProfiler,",",fixed=T)[[1]]

	# For all genes
	for (CluterProfiler_index in all_genes_of_annotation)
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
		gene_annotation<-data.frame(gene_id=gene_id,gene=gene,log2change=log2change,Category=Category,Pvalue=Pvalue,FDR=FDR,ENTREZID=ENTREZID,Symbol=Symbol,Description=Description,genecards_Category=genecards_Category,UniProt_ID=UniProt_ID,GIFtS=GIFtS,GC_id=GC_id,GeneCards_Summary=GeneCards_Summary,Stage=Stage,CluterProfiler=CluterProfiler_index )

		# df_all_annotation
		df_all_annotation<-rbind(df_all_annotation,gene_annotation)
	}	
}
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
n_top<-2

# Stage I
Stage="Stage I"

# Store annotation of stage I
df_all_annotation_per_stage<-df_all_annotation[which(df_all_annotation$Stage==Stage),]
######################################################################################################################
# A table for the count of go terms
table_GO<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="GO","Description"])
table_KEGG<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="KEGG","Description"])
table_REACTOME<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="Reactome","Description"])

selection_GO_Stage_I          <-names(tail(sort(table_GO),n=n_top))
selection_KEGG_Stage_I        <-names(tail(sort(table_KEGG),n=n_top))
selection_REACTOM_Stage_I     <-names(tail(sort(table_REACTOME),n=n_top))
######################################################################################################################
# Stage II
Stage="Stage II"

# Store annotation of stage II
df_all_annotation_per_stage<-df_all_annotation[which(df_all_annotation$Stage==Stage),]
######################################################################################################################
# A table for the count of go terms
table_GO<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="GO","Description"])
table_KEGG<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="KEGG","Description"])
table_REACTOME<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="Reactome","Description"])

selection_GO_Stage_II          <-names(tail(sort(table_GO),n=n_top))
selection_KEGG_Stage_II        <-names(tail(sort(table_KEGG),n=n_top))
selection_REACTOM_Stage_II     <-names(tail(sort(table_REACTOME),n=n_top))
######################################################################################################################
# Stage III
Stage="Stage III"

# Store annotation of stage III
df_all_annotation_per_stage<-df_all_annotation[which(df_all_annotation$Stage==Stage),]
######################################################################################################################
# A table for the count of go terms
table_GO<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="GO","Description"])
table_KEGG<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="KEGG","Description"])
table_REACTOME<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="Reactome","Description"])

selection_GO_Stage_III          <-names(tail(sort(table_GO),n=n_top))
selection_KEGG_Stage_III        <-names(tail(sort(table_KEGG),n=n_top))
selection_REACTOM_Stage_III     <-names(tail(sort(table_REACTOME),n=n_top))
######################################################################################################################

######################################################################################################################
# Store annotation of stages
df_all_annotation_per_stage<-df_all_annotation

# interactome_annotation
interactome_annotation_stage<-interactome_annotation

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
gene_Stage_I                  <-unique(df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage I"),"SYMBOL"])
gene_Stage_II                 <-unique(df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage II"),"SYMBOL"])
gene_Stage_III                 <-unique(df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage III"),"SYMBOL"])

selection_GO          <-unique(c(selection_GO_Stage_I,selection_GO_Stage_II,selection_GO_Stage_III))
selection_KEGG        <-unique(c(selection_KEGG_Stage_I,selection_KEGG_Stage_II,selection_KEGG_Stage_III))
selection_REACTOM     <-unique(c(selection_REACTOM_Stage_I,selection_REACTOM_Stage_II,selection_REACTOM_Stage_III))

df_all_annotation_selected_pathways<-rbind(df_all_annotation_per_stage[which(df_all_annotation_per_stage$Description %in% selection_GO),],
df_all_annotation_per_stage[which(df_all_annotation_per_stage$Description %in% selection_KEGG),],
df_all_annotation_per_stage[which(df_all_annotation_per_stage$Description %in% selection_REACTOM),],
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
V(graph_all_stages)$color                                                                         <- "black"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_GO)]$color       <- "gainsboro"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_KEGG)]$color     <- "lightslategray"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_REACTOME)]$color <- "darkslategray"

V(graph_all_stages)[which(names(V(graph_all_stages)) %in% gene_Stage_I)]$color       <- "#E69F00"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% gene_Stage_II)]$color      <- "#009E73"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% gene_Stage_III)]$color     <- "#D81B60"


# Vertice colours of genes
V(graph_all_stages)$shape                                                                           <-"circle"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_GO)]$shape        <- "square"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_KEGG)]$shape      <- "square"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_REACTOME)]$shape  <- "square"

# Set size of the node according to the dregree
V(graph_all_stages)$size                                                                          <- 5
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_GO)]$size        <- 7
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_KEGG)]$size      <- 7
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% annotation_description_REACTOME)]$size  <- 7

# Vertice colours of genes
E(graph_all_stages)$color                                                                         <- "black"

# Set ronames
E(graph_all_stages)[which(df_all_eges$names %in% selected_interactome)]$color                     <- "lightgrey"

# FindClusters_resolution
png(filename=paste(output_folder,"Plot_Stage_all.png",sep=""), width = 30, height = 30, res=600, units = "cm")
	#plot(graph_all_stages, layout=layout_nicely, vertex.label=NA) # Stage I
	#plot(graph_all_stages, layout=  layout_with_kk, vertex.label=NA) # Stage II
	plot(graph_all_stages, layout=   layout_nicely, vertex.label=NA) # Stage II
dev.off()

# FindClusters_resolution
png(filename=paste(output_folder,"Plot_Stage_label.png",sep=""), width = 30, height = 30, res=600, units = "cm")
	#plot(graph_all_stages, layout=layout_nicely) # Stage I
	#plot(graph_all_stages, layout=  layout_with_kk) # Stage II
	plot(graph_all_stages, layout=   layout_nicely) # Stage II
dev.off()


# Save file 
write.xlsx(x=df_all_annotation_selected_pathways,file=paste(output_dir,"df_all_annotation",".xlsx",sep=""), sheet="Stage I", append=FALSE)


# Set legend
png(filename=paste(output_folder,"legend.png",sep=""), width = 5, height = 5, res=600, units = "cm")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c('GO', 'KEGG', 'REACTOME','Gene'), pch=16, pt.cex=3, cex=1.5, bty='n',col = c('#ff6347', '#ffd700', '#7f7f7f', '#0072b2'))
mtext("Legend", at=0.2, cex=2)
dev.off()
