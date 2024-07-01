######################################################################################################################
# Specify sheet by its name
# Read the table with manual and automatic annotation
annotation_stages_all <- data.frame(read_excel("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/unique_genes_annotation_clusterProfiler_Genecards.xlsx"))
######################################################################################################################
# Specify sheet by its name
# split annotation per stage
annotation_stage_I <- data.frame(read_excel("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/unique_genes_annotation.xlsx", sheet = "Stage I"))
annotation_stage_II <- data.frame(read_excel("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/unique_genes_annotation.xlsx", sheet = "Stage II"))
annotation_stage_III <- data.frame(read_excel("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/unique_genes_annotation.xlsx", sheet = "Stage III"))
######################################################################################################################
# A data.frame to store all results
# Save results per annotation term, below only annotation of clusterProfiler
df_all_annotation<-data.frame(gene_id=c(),gene=c(),log2change=c(),Category=c(),Pvalue=c(),FDR=c(),ENTREZID=c(),Symbol=c(),Description=c(),genecards_Category=c(),UniProt_ID=c(),GIFtS=c(),GC_id=c(),GeneCards_Summary=c(),Stage=c(),Layer=c(),CluterProfiler=c() )

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
	GeneCards	     <- annotation_stages_all[annotation,"GeneCards"]
	Layer       	     <- ""

	# All genes of stage
	all_CluterProfiler<-strsplit(x=CluterProfiler,",",fixed=T)[[1]]
	all_GeneCards<-strsplit(x=GeneCards,",",fixed=T)[[1]]
	print(all_CluterProfiler)
	
	# For all genes
	for (CluterProfiler_index in all_CluterProfiler)
	{
		# Set layer
		Layer<-strsplit(x=all_CluterProfiler,":",fixed=T)[[1]][[1]]
		
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
		gene_annotation<-data.frame(gene_id=gene_id,gene=gene,log2change=log2change,Category=Category,Pvalue=Pvalue,FDR=FDR,ENTREZID=ENTREZID,Symbol=Symbol,Description=Description,genecards_Category=genecards_Category,UniProt_ID=UniProt_ID,GIFtS=GIFtS,GC_id=GC_id,GeneCards_Summary=GeneCards_Summary,Stage=Stage,Layer=Layer,CluterProfiler=str_trim(CluterProfiler_index, side="left") )

		# df_all_annotation
		df_all_annotation<-rbind(df_all_annotation,gene_annotation)
	}	
}
######################################################################################################################
# The ten most abundant annotation terms in number of genes were selected for further inspection (see Figure Annotaton): Reactome:Intracellular signaling by second messengers (14), Reactome:PIP3 activates AKT signaling (14), Reactome:Regulation of expression of SLITs and ROBOs  (14), Reactome:Signaling by ROBO receptors (15), GO:mitochondrial matrix (16), GO:mitochondrial protein-containing complex (16), KEGG:Alzheimer disease (16), Reactome:Translation  (16), GO:mitochondrial inner membrane (17), Reactome:SARS-CoV Infections (18)
######################################################################################################################
# The ten most abundant annotation terms in number of genes were selected for further inspection (see Figure Annotaton): Reactome:Intracellular signaling by second messengers (14), Reactome:PIP3 activates AKT signaling (14), Reactome:Regulation of expression of SLITs and ROBOs  (14), Reactome:Signaling by ROBO receptors (15), GO:mitochondrial matrix (16), GO:mitochondrial protein-containing complex (16), KEGG:Alzheimer disease (16), Reactome:Translation  (16), GO:mitochondrial inner membrane (17), Reactome:SARS-CoV Infections (18)
# Load interactome from file
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadInteractomeUniqueGenes.R")

# gene annotation
interactome_annotation<-data.frame(gene_id=rownames(interactome_all_stage),gene=rownames(interactome_all_stage),log2change=rownames(interactome_all_stage),Category=0,Pvalue=0,FDR=0,ENTREZID=0,Symbol=interactome_all_stage$Gene1,Description=0,genecards_Category=0,UniProt_ID=0,GIFtS=0,GC_id=0,GeneCards_Summary=0,Stage=interactome_all_stage$Stage,Layer="Interactome",CluterProfiler=interactome_all_stage$Gene2 )
coexpression_annotation<-data.frame(gene_id=rownames(net_stage_all_correlation_network),gene=rownames(net_stage_all_correlation_network),log2change=rownames(net_stage_all_correlation_network),Category=0,Pvalue=0,FDR=0,ENTREZID=0,Symbol=net_stage_all_correlation_network$Gene1,Description=0,genecards_Category=0,UniProt_ID=0,GIFtS=0,GC_id=0,GeneCards_Summary=0,Stage=net_stage_all_correlation_network$Stage,Layer="Coexpression",CluterProfiler=net_stage_all_correlation_network$Gene2 )
######################################################################################################################
selected_interactome<-paste(interactome_all_stage$Gene1,interactome_all_stage$Gene2, sep="-")
selected_coexpression<-paste(net_stage_all_correlation_network$Gene1,net_stage_all_correlation_network$Gene2, sep="-")

# Next,  ten most abundant annotation terms in number of genes per stage additionally inspected.
######################################################################################################################
# Store annotation of stages
df_all_annotation_per_stage<-df_all_annotation

# interactome_annotation
interactome_annotation_stage<-interactome_annotation
######################################################################################################################
# A table for the count of go terms
table_Stage_I   <-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Stage=="Stage I","CluterProfiler"])
table_Stage_II  <-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Stage=="Stage II","CluterProfiler"])
table_Stage_III <-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Stage=="Stage III","CluterProfiler"])

stage_I_anotation<-df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage I"),]
stage_II_anotation<-df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage II"),]
stage_III_anotation<-df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage III"),]

selection<-c(names(table_Stage_I),names(table_Stage_II),names(table_Stage_III))

# df_count_terms
df_count_terms<-data.frame(Term=unique(selection),Stage_I=0,Stage_II=0,Stage_III=0)

# Set rownames
rownames(df_count_terms)<-df_count_terms$Term

# for each term, count number of genes
for (term in df_count_terms$Term)
{	
	# Save counts
	df_count_terms[term,"Stage_I"]   <-length(stage_I_anotation[stage_I_anotation$CluterProfiler==term,"Symbol"])
	df_count_terms[term,"Stage_II"]  <-length(stage_II_anotation[stage_II_anotation$CluterProfiler==term,"Symbol"])
	df_count_terms[term,"Stage_III"] <-length(stage_III_anotation[stage_III_anotation$CluterProfiler==term,"Symbol"])
}
df_count_terms<-df_count_terms[rownames(df_count_terms)!="-",]
###########################################################################################################################
# Save file 
write.xlsx(x=df_count_terms,file=paste(output_dir,"unique_genes_annotation_count",".xlsx",sep=""), sheet="genes and terms per stage", append=FALSE)


# Selection terms to be analyses
# I stopped here

# Boolean count per terms
boolean_counts_terms<-df_count_terms[,2:4]>0
######################################################################################################################
#table_Stage_I<-table_Stage_I[which(table_Stage_I>1)]
#table_Stage_II<-table_Stage_II[which(table_Stage_II>1)]
#table_Stage_III<-table_Stage_III[which(table_Stage_III>1)]
####################################################################################################################
# Store information for each gene
df_count_terms_selected<-df_count_terms

# Genes for each stage
df_count_terms_selected$Genes_Stage_I<-""
df_count_terms_selected$Genes_Stage_II<-""
df_count_terms_selected$Genes_Stage_III<-""


# for each term, take the annotation about the genes
for (term in rownames(df_count_terms_selected))
{
	# df_all_annotation_stage
	df_all_annotation_stage_I  <-df_all_annotation[df_all_annotation$Stage=="Stage I",]
	df_all_annotation_stage_II <-df_all_annotation[df_all_annotation$Stage=="Stage II",]
	df_all_annotation_stage_III<-df_all_annotation[df_all_annotation$Stage=="Stage III",]

	df_count_terms_selected[term,"Genes_Stage_I"]<-paste(df_all_annotation_stage_I[df_all_annotation_stage_I$CluterProfiler==term,"Symbol"],collapse=", ")
	df_count_terms_selected[term,"Genes_Stage_II"]<-paste(df_all_annotation_stage_II[df_all_annotation_stage_II$CluterProfiler==term,"Symbol"],collapse=", ")
	df_count_terms_selected[term,"Genes_Stage_III"]<-paste(df_all_annotation_stage_III[df_all_annotation_stage_III$CluterProfiler==term,"Symbol"],collapse=", ")
}
# Save file 
write.xlsx(x=df_count_terms_selected,file=paste(output_dir,"unique_genes_annotation_count",".xlsx",sep=""), sheet="all pathways", append=TRUE)
####################################################################################################################
# To do : 1st July 2024
# For each of the term, take the stage-specific genes associated to it per stage.
df_count_terms_selected$Stage_I_norm   <- df_count_terms_selected$Stage_I/dim(annotation_stage_I)[1]*100
df_count_terms_selected$Stage_II_norm  <- df_count_terms_selected$Stage_II/dim(annotation_stage_II)[1]*100     
df_count_terms_selected$Stage_III_norm <- df_count_terms_selected$Stage_III/dim(annotation_stage_III)[1]*100

####################################################################################################################
# Here I will split the count by category  GO, Reactome, KEGG
# To do : 1st July 2024
# Split into "GO", "Reactome" and "KEGG" based of the string of the term name containts these words
df_count_terms_selected_GO       <- data.frame(Term=c(),Stage_I=c(),Stage_II=c(),Stage_III=c(),Genes_Stage_I=c(), Genes_Stage_II=c(), Genes_Stage_III=c())
df_count_terms_selected_Reactome <- data.frame(Term=c(),Stage_I=c(),Stage_II=c(),Stage_III=c(),Genes_Stage_I=c(), Genes_Stage_II=c(), Genes_Stage_III=c())
df_count_terms_selected_KEGG     <- data.frame(Term=c(),Stage_I=c(),Stage_II=c(),Stage_III=c(),Genes_Stage_I=c(), Genes_Stage_II=c(), Genes_Stage_III=c())

# For each row of the the counts tanble, check if belong to "GO", "Reactome" and "KEGG"
for (selected_term in rownames(df_count_terms_selected))
{
	# If contains the word GO
	if(grepl("GO:", selected_term, fixed = TRUE))
	{	
		# Add to the data.frame of GO
		df_count_terms_selected_GO<-rbind(df_count_terms_selected_GO, df_count_terms_selected[selected_term,])
	}
	# If contains the word GO
	if(grepl("KEGG:", selected_term, fixed = TRUE))
	{	
		# Add to the data.frame of GO
		df_count_terms_selected_KEGG<-rbind(df_count_terms_selected_KEGG, df_count_terms_selected[selected_term,])
	}	
	# If contains the word GO
	if(grepl("Reactome:", selected_term, fixed = TRUE))
	{	
		# Add to the data.frame of GO
		df_count_terms_selected_Reactome<-rbind(df_count_terms_selected_Reactome, df_count_terms_selected[selected_term,])
	}		
}
####################################################################################################################
write.xlsx(x=df_count_terms_selected_GO,file=paste(output_dir,"unique_genes_annotation_count",".xlsx",sep=""), sheet="GO terms", append=TRUE)
write.xlsx(x=df_count_terms_selected_KEGG,file=paste(output_dir,"unique_genes_annotation_count",".xlsx",sep=""), sheet="KEGG terms", append=TRUE)
write.xlsx(x=df_count_terms_selected_Reactome,file=paste(output_dir,"unique_genes_annotation_count",".xlsx",sep=""), sheet="Reactome terms", append=TRUE)
####################################################################################################################
# Remove empty line
# Select top 10 terms                                                                                                                                                            #
selection_GO       <-unique(c(names(tail(sort(table_GO_Stage_I),n=10)),names(tail(sort(table_GO_Stage_II),n=10)),names(tail(sort(table_GO_Stage_III),n=10))))                    #
selection_Reactome <-unique(c(names(tail(sort(table_Reactome_Stage_I),n=10)),names(tail(sort(table_Reactome_Stage_II),n=10)),names(tail(sort(table_Reactome_Stage_III),n=10))))  #
selection_KEGG     <-unique(c(names(tail(sort(table_KEGG_Stage_I),n=10)),names(tail(sort(table_KEGG_Stage_II),n=10)),names(tail(sort(table_KEGG_Stage_III),n=10))))              #

# Store information for each gene
matrix_count_terms_selected_GO        <-df_count_terms_selected_GO[selection_GO,]
matrix_count_terms_selected_KEGG      <-df_count_terms_selected_KEGG[selection_KEGG,]
matrix_count_terms_selected_Reactome  <-df_count_terms_selected_Reactome[selection_Reactome,]

go_order  <-hcluster(matrix_count_terms_selected_GO[,c("Stage_I_norm","Stage_II_norm","Stage_III_norm")],link = "ave")$labels[hcluster(matrix_count_terms_selected_GO[,c("Stage_I_norm","Stage_II_norm","Stage_III_norm")],link = "ave")$order]
kegg_order<-hcluster(matrix_count_terms_selected_KEGG[,c("Stage_I_norm","Stage_II_norm","Stage_III_norm")],link = "ave")$labels[hcluster(matrix_count_terms_selected_KEGG[,c("Stage_I_norm","Stage_II_norm","Stage_III_norm")],link = "ave")$order]
reactome_order<-hcluster(matrix_count_terms_selected_Reactome[,c("Stage_I_norm","Stage_II_norm","Stage_III_norm")],link = "ave")$labels[hcluster(matrix_count_terms_selected_Reactome[,c("Stage_I_norm","Stage_II_norm","Stage_III_norm")],link = "ave")$order]

matrix_count_terms_selected_GO  <-matrix_count_terms_selected_GO[go_order,]
matrix_count_terms_selected_KEGG<-matrix_count_terms_selected_KEGG[kegg_order,]
matrix_count_terms_selected_Reactome<-matrix_count_terms_selected_Reactome[reactome_order,]
####################################################################################################################
# Save file 
write.xlsx(x=matrix_count_terms_selected_GO,file=paste(output_dir,"unique_genes_annotation_count",".xlsx",sep=""), sheet="selected GO", append=TRUE)
write.xlsx(x=matrix_count_terms_selected_KEGG,file=paste(output_dir,"unique_genes_annotation_count",".xlsx",sep=""), sheet="selected KEGG", append=TRUE)
write.xlsx(x=matrix_count_terms_selected_Reactome,file=paste(output_dir,"unique_genes_annotation_count",".xlsx",sep=""), sheet="selected Reactome", append=TRUE)
####################################################################################################################









####################################################################################################################
# Concatenate table
df_all_annotation_selected_pathways<-rbind(df_all_annotation[which(df_all_annotation$CluterProfiler %in% selection_GO),],
interactome_annotation_stage,coexpression_annotation)
####################################################################################################################
# All stages
graph_all_stages <- graph_from_data_frame(d=unique(df_all_annotation_selected_pathways[,c("Symbol","CluterProfiler")]), vertices=unique(c(df_all_annotation_selected_pathways$Symbol,df_all_annotation_selected_pathways$CluterProfiler)), directed=F)  

# df_all_eges
df_all_eges<-data.frame(get.edgelist(graph_all_stages))

# Set ronames
df_all_eges$names<-paste(df_all_eges$X1,df_all_eges$X2,sep="-")
####################################################################################################################
# Vertice colours of genes
V(graph_all_stages)$color                                                                                                  <- "black"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% df_all_annotation_selected_pathways$CluterProfiler)]$color       <- "black"

stage_I   <-names(V(graph_all_stages)[which(names(V(graph_all_stages)) %in% ids_stage_I$SYMBOL)])
stage_II  <-names(V(graph_all_stages)[which(names(V(graph_all_stages)) %in% ids_stage_II$SYMBOL)])
stage_III <-names(V(graph_all_stages)[which(names(V(graph_all_stages)) %in% ids_stage_III$SYMBOL)])

names(V(graph_all_stages)) %in% c(stage_I,stage_II,stage_III)
length(names(V(graph_all_stages)) %in% c(stage_I,stage_II,stage_III))
names(V(graph_all_stages))[!names(V(graph_all_stages)) %in% c(stage_I,stage_II,stage_III)]


V(graph_all_stages)[which(names(V(graph_all_stages)) %in% ids_stage_I$SYMBOL)]$color       <- "#E69F00"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% ids_stage_II$SYMBOL)]$color      <- "#009E73"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% ids_stage_III$SYMBOL)]$color     <- "#D81B60"

# Vertice colours of genes
V(graph_all_stages)$shape                                                                                                                                                                                      <-"circle"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in%  df_all_annotation_selected_pathways[df_all_annotation_selected_pathways$Layer %in% c("GO","KEGG","Reactome"),"CluterProfiler"])]$shape              <- "square"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in%  df_all_annotation_selected_pathways[df_all_annotation_selected_pathways$Layer %in% c("Ontology","Disease","Pathway"),"CluterProfiler"])]$shape      <- "square"

# Set size of the node according to the dregree
V(graph_all_stages)$size                                                                                                                                                                                       <- 5
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% df_all_annotation_selected_pathways[df_all_annotation_selected_pathways$Layer %in% c("Ontology","Disease","Pathway"),"CluterProfiler"])]$size        <- 7
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% df_all_annotation_selected_pathways[df_all_annotation_selected_pathways$Layer %in% c("GO","KEGG","Reactome"),"CluterProfiler"])]$size                <- 7

# Vertice colours of genes
E(graph_all_stages)$color                                                                         <- "lightgrey"

# Set ronames
E(graph_all_stages)[which(df_all_eges$names %in% selected_interactome)]$color                     <- "black"
E(graph_all_stages)[which(df_all_eges$names %in% selected_coexpression)]$color                    <- "darkblue"
matrix_count_terms_selected_GO        <-df_count_terms_selected_GO[selection_GO,]
matrix_count_terms_selected_KEGG      <-df_count_terms_selected_KEGG[selection_KEGG,]
matrix_count_terms_selected_Reactome  <-df_count_terms_selected_Reactome[selection_Reactome,]

# Set size of the node according to the dregree
V(graph_all_stages)$label                                                                                                                                                                                      <- ""
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% df_all_annotation_selected_pathways[df_all_annotation_selected_pathways$Layer %in% c("Ontology","Disease","Pathway"),"CluterProfiler"])]$label        <- names(V(graph_all_stages)[which(names(V(graph_all_stages)) %in% df_all_annotation_selected_pathways[df_all_annotation_selected_pathways$Layer %in% c("Ontology","Disease","Pathway"),"CluterProfiler"])])
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% df_all_annotation_selected_pathways[df_all_annotation_selected_pathways$Layer %in% c("GO","KEGG","Reactome"),"CluterProfiler"])]$label                <- names(V(graph_all_stages)[which(names(V(graph_all_stages)) %in% df_all_annotation_selected_pathways[df_all_annotation_selected_pathways$Layer %in% c("GO","KEGG","Reactome"),"CluterProfiler"])])


go_order  <-hcluster(matrix_count_terms_selected_GO[,c("Stage_I_norm","Stage_II_norm","Stage_III_norm")],link = "ave")$labels[hcluster(matrix_count_terms_selected_GO[,c("Stage_I_norm","Stage_II_norm","Stage_III_norm")],link = "ave")$order]
kegg_order<-hcluster(matrix_count_terms_selected_KEGG[,c("Stage_I_norm","Stage_II_norm","Stage_III_norm")],link = "ave")$labels[hcluster(matrix_count_terms_selected_KEGG[,c("Stage_I_norm","Stage_II_norm","Stage_III_norm")],link = "ave")$order]
reactome_order<-hcluster(matrix_count_terms_selected_Reactome[,c("Stage_I_norm","Stage_II_norm","Stage_III_norm")],link = "ave")$labels[hcluster(matrix_count_terms_selected_Reactome[,c("Stage_I_norm","Stage_II_norm","Stage_III_norm")],link = "ave")$order]

# FindClusters_resolution
png(filename=paste(output_folder,"Plot_Stage_all_per_Stagge.png",sep=""), width = 30, height = 30, res=600, units = "cm")
	#plot(graph_all_stages, layout=layout_nicely, vertex.label=NA) # Stage I
	#plot(graph_all_stages, layout=  layout_with_kk, vertex.label=NA) # Stage II
	#plot(graph_all_stages, layout=   layout_nicely, vertex.label=V(graph_all_stages)$labels) # Stage II
	plot(graph_all_stages, layout=   layout_nicely,vertex.label=V(graph_all_stages)$label, edge.colour=NA ) # Stage II
dev.off()






###########################################################################
# FindClusters_resolution
png(filename=paste(output_folder,"Plot_Stage_label.png",sep=""), width = 30, height = 30, res=600, units = "cm")
	#plot(graph_all_stages, layout=layout_nicely) # Stage I
	#plot(graph_all_stages, layout=  layout_with_kk) # Stage II
	plot(graph_all_stages, layout=   layout_nicely) # Stage II
dev.off()


# Save file 
write.xlsx(x=df_all_annotation_selected_pathways,file=paste(output_dir,"unique_genes_annotation_clusterProfiler",".xlsx",sep=""), sheet="ten most abundant per stage", append=TRUE)

# Set legend
png(filename=paste(output_folder,"legend.png",sep=""), width = 5, height = 5, res=600, units = "cm")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c('GO', 'KEGG', 'REACTOME','Gene'), pch=16, pt.cex=3, cex=1.5, bty='n',col = c('#ff6347', '#ffd700', '#7f7f7f', '#0072b2'))
mtext("Legend", at=0.2, cex=2)
###########################################################################
