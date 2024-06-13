annotation_stages_all
######################################################################################################################
# A data.frame to store all results
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
	Layer       	     <- ""

	# All genes of stage
	all_genes_of_annotation<-strsplit(x=CluterProfiler,",",fixed=T)[[1]]

	# For all genes
	for (CluterProfiler_index in all_genes_of_annotation)
	{
		# Set layer
		Layer<-strsplit(x=CluterProfiler_index,":",fixed=T)[[1]][[1]]
		
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
		gene_annotation<-data.frame(gene_id=gene_id,gene=gene,log2change=log2change,Category=Category,Pvalue=Pvalue,FDR=FDR,ENTREZID=ENTREZID,Symbol=Symbol,Description=Description,genecards_Category=genecards_Category,UniProt_ID=UniProt_ID,GIFtS=GIFtS,GC_id=GC_id,GeneCards_Summary=GeneCards_Summary,Stage=Stage,Layer=Layer,CluterProfiler=CluterProfiler_index )

		# df_all_annotation
		df_all_annotation<-rbind(df_all_annotation,gene_annotation)
	}	
}












######################################################################################################################
# The ten most abundant annotation terms in number of genes were selected for further inspection (see Figure Annotaton): Reactome:Intracellular signaling by second messengers (14), Reactome:PIP3 activates AKT signaling (14), Reactome:Regulation of expression of SLITs and ROBOs  (14), Reactome:Signaling by ROBO receptors (15), GO:mitochondrial matrix (16), GO:mitochondrial protein-containing complex (16), KEGG:Alzheimer disease (16), Reactome:Translation  (16), GO:mitochondrial inner membrane (17), Reactome:SARS-CoV Infections (18)
######################################################################################################################
# The ten most abundant annotation terms in number of genes were selected for further inspection (see Figure Annotaton): Reactome:Intracellular signaling by second messengers (14), Reactome:PIP3 activates AKT signaling (14), Reactome:Regulation of expression of SLITs and ROBOs  (14), Reactome:Signaling by ROBO receptors (15), GO:mitochondrial matrix (16), GO:mitochondrial protein-containing complex (16), KEGG:Alzheimer disease (16), Reactome:Translation  (16), GO:mitochondrial inner membrane (17), Reactome:SARS-CoV Infections (18)
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
interactome_annotation<-data.frame(gene_id=rownames(interactome_all_stage),gene=rownames(interactome_all_stage),log2change=rownames(interactome_all_stage),Category=0,Pvalue=0,FDR=0,ENTREZID=0,Symbol=interactome_all_stage$Gene1,Description=0,genecards_Category=0,UniProt_ID=0,GIFtS=0,GC_id=0,GeneCards_Summary=0,Stage=interactome_all_stage$Stage,Layer="Interactome",CluterProfiler=interactome_all_stage$Gene2 )
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
table_GO<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="GO","CluterProfiler"])
table_KEGG<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="KEGG","CluterProfiler"])
table_REACTOME<-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Layer=="Reactome","CluterProfiler"])

all_selection<-c(table_GO,table_KEGG,table_REACTOME)

selection_all         <-names(tail(sort(all_selection),n=10))

# Count per selected term
data.frame(n=tail(sort(all_selection),n=10))

df_all_annotation_selected_pathways<-rbind(df_all_annotation_per_stage[which(df_all_annotation_per_stage$CluterProfiler %in% selection_all),],
interactome_annotation_stage)
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

V(graph_all_stages)[which(names(V(graph_all_stages)) %in% ids_stage_I$SYMBOL)]$color       <- "#E69F00"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% ids_stage_II$SYMBOL)]$color      <- "#009E73"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% ids_stage_III$SYMBOL)]$color     <- "#D81B60"

# Vertice colours of genes
V(graph_all_stages)$shape                                                                                                                                                                                      <-"circle"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in%  df_all_annotation_selected_pathways[df_all_annotation_selected_pathways$Layer %in% c("GO","KEGG","Reactome"),"CluterProfiler"])]$shape              <- "square"

# Set size of the node according to the dregree
V(graph_all_stages)$size                                                                                                                                                                               <- 5
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% df_all_annotation_selected_pathways[df_all_annotation_selected_pathways$Layer %in% c("GO","KEGG","Reactome"),"CluterProfiler"])]$size        <- 7

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
write.xlsx(x=df_all_annotation_selected_pathways,file=paste(output_dir,"unique_genes_annotation_clusterProfiler",".xlsx",sep=""), sheet="ten most abundant annotation", append=TRUE)

# Set legend
png(filename=paste(output_folder,"legend.png",sep=""), width = 5, height = 5, res=600, units = "cm")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c('GO', 'KEGG', 'REACTOME','Gene'), pch=16, pt.cex=3, cex=1.5, bty='n',col = c('#ff6347', '#ffd700', '#7f7f7f', '#0072b2'))
mtext("Legend", at=0.2, cex=2)
dev.off()
















######################################################################################################################
# Next,  ten most abundant annotation terms in number of genes per stage additionally inspected.
######################################################################################################################
# Store annotation of stages
df_all_annotation_per_stage<-df_all_annotation

# interactome_annotation
interactome_annotation_stage<-interactome_annotation

# interactome_annotation_stage
selected_interactome<-c(paste(interactome_annotation_stage$Description,interactome_annotation_stage$SYMBOL,sep="-"),
paste(interactome_annotation_stage$SYMBOL,interactome_annotation_stage$Description,sep="-"))
######################################################################################################################
# Se
# A table for the count of go terms
table_Stage_I   <-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Stage=="Stage I","CluterProfiler"])
table_Stage_II  <-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Stage=="Stage II","CluterProfiler"])
table_Stage_III <-table(df_all_annotation_per_stage[df_all_annotation_per_stage$Stage=="Stage III","CluterProfiler"])

table_Stage_I<-table_Stage_I[which(table_Stage_I>1)]
table_Stage_II<-table_Stage_II[which(table_Stage_II>1)]
table_Stage_II<-table_Stage_III[which(table_Stage_III>1)]

# select terms
selection_all <-unique(c(names(tail(sort(table_Stage_I),n=20)),names(tail(sort(table_Stage_II),n=20)),names(tail(sort(table_Stage_III),n=20))))

stage_I_anotation<-df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage I"),]
stage_II_anotation<-df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage II"),]
stage_III_anotation<-df_all_annotation_per_stage[which(df_all_annotation_per_stage$Stage=="Stage III"),]

# df_count_terms
df_count_terms<-data.frame(Term=selection_all,Stage_I=0,Stage_II=0,Stage_III=0)

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
###########################################################################################################################
df_count_terms$Stage_I_Genes<-""
df_count_terms$Stage_II_Genes<-""
df_count_terms$Stage_III_Genes<-""

# A table for the count of go terms
table_Stage_I     <-df_all_annotation_per_stage[df_all_annotation_per_stage$Stage=="Stage I",]
table_Stage_II    <-df_all_annotation_per_stage[df_all_annotation_per_stage$Stage=="Stage II",]
table_Stage_III   <-df_all_annotation_per_stage[df_all_annotation_per_stage$Stage=="Stage III",]


# For each term 
for (term in rownames(df_count_terms))
{
	df_count_terms[term,"Stage_I_Genes"]  <-paste(stage_I_anotation[stage_I_anotation$CluterProfiler==term,"Symbol"],collapse=",")
	df_count_terms[term,"Stage_II_Genes"] <-paste(stage_II_anotation[stage_II_anotation$CluterProfiler==term,"Symbol"],collapse=",")
	df_count_terms[term,"Stage_III_Genes"]<-paste(stage_III_anotation[stage_III_anotation$CluterProfiler==term,"Symbol"],collapse=",")
}

# Save file 
write.xlsx(x=df_count_terms,file=paste(output_dir,"unique_genes_annotation_clusterProfiler",".xlsx",sep=""), sheet="genes and terms per stage", append=FALSE)

# Boolean count per terms
boolean_counts_terms<-df_count_terms[,2:4]>1
###########################################################################################################################
# Question 1: terms in all stages
paste(names(which(boolean_counts_terms[,c("Stage_I")] & boolean_counts_terms[,c("Stage_II")] & boolean_counts_terms[,c("Stage_III")])),collapse=", ")
# Answer : Reactome:Asparagine N-linked glycosylation, Reactome:Neutrophil degranulation
# Genes in these terms
names(which(boolean_counts_terms[,c("Stage_I")] & boolean_counts_terms[,c("Stage_II")] & boolean_counts_terms[,c("Stage_III")]))
paste(df_all_annotation_per_stage[df_all_annotation_per_stage$CluterProfiler %in% names(which(boolean_counts_terms[,c("Stage_I")] & boolean_counts_terms[,c("Stage_II")] & boolean_counts_terms[,c("Stage_III")])),"Symbol"], collapse=", ")
# paste(df_all_annotation_per_stage[df_all_annotation_per_stage$CluterProfiler %in% names(which(boolean_counts_terms[,c("Stage_I")] & boolean_counts_terms[,c("Stage_II")] & boolean_counts_terms[,c("Stage_III")])),"Symbol"], collapse=", ")
######################################################################################################################

df_all_annotation_selected_pathways<-rbind(df_all_annotation_per_stage[which(df_all_annotation_per_stage$CluterProfiler %in% selection_all),],
interactome_annotation_stage)
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

V(graph_all_stages)[which(names(V(graph_all_stages)) %in% ids_stage_I$SYMBOL)]$color       <- "#E69F00"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% ids_stage_II$SYMBOL)]$color      <- "#009E73"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% ids_stage_III$SYMBOL)]$color     <- "#D81B60"

# Vertice colours of genes
V(graph_all_stages)$shape                                                                                                                                                                                      <-"circle"
V(graph_all_stages)[which(names(V(graph_all_stages)) %in%  df_all_annotation_selected_pathways[df_all_annotation_selected_pathways$Layer %in% c("GO","KEGG","Reactome"),"CluterProfiler"])]$shape              <- "square"

# Set size of the node according to the dregree
V(graph_all_stages)$size                                                                                                                                                                               <- 5
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% df_all_annotation_selected_pathways[df_all_annotation_selected_pathways$Layer %in% c("GO","KEGG","Reactome"),"CluterProfiler"])]$size        <- 7

# Vertice colours of genes
E(graph_all_stages)$color                                                                         <- "black"

# Set ronames
E(graph_all_stages)[which(df_all_eges$names %in% selected_interactome)]$color                     <- "lightgrey"

# FindClusters_resolution
png(filename=paste(output_folder,"Plot_Stage_all_per_Stagge.png",sep=""), width = 30, height = 30, res=600, units = "cm")
	#plot(graph_all_stages, layout=layout_nicely, vertex.label=NA) # Stage I
	#plot(graph_all_stages, layout=  layout_with_kk, vertex.label=NA) # Stage II
	plot(graph_all_stages, layout=   layout_nicely, vertex.label=NA) # Stage II
dev.off()

###########################################################################3
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

