# Specify sheet by its name
annotation_stages_all <- data.frame(read_excel("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/unique_genes_annotation_clusterProfiler_Genecards.xlsx"))
######################################################################################################################
# Specify sheet by its name
annotation_stage_I <- data.frame(read_excel("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/unique_genes_annotation.xlsx", sheet = "Stage I"))
annotation_stage_II <- data.frame(read_excel("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/unique_genes_annotation.xlsx", sheet = "Stage II"))
annotation_stage_III <- data.frame(read_excel("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/unique_genes_annotation.xlsx", sheet = "Stage III"))
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
	GeneCards	     <- annotation_stages_all[annotation,"GeneCards"]
	Layer       	     <- ""

	# All genes of stage
	all_CluterProfiler<-strsplit(x=CluterProfiler,",",fixed=T)[[1]]
	all_GeneCards<-strsplit(x=GeneCards,",",fixed=T)[[1]]
	
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
	# For all genes
	for (Genecards_index in all_GeneCards)
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
		gene_annotation<-data.frame(gene_id=gene_id,gene=gene,log2change=log2change,Category=Category,Pvalue=Pvalue,FDR=FDR,ENTREZID=ENTREZID,Symbol=Symbol,Description=Description,genecards_Category=genecards_Category,UniProt_ID=UniProt_ID,GIFtS=GIFtS,GC_id=GC_id,GeneCards_Summary=GeneCards_Summary,Stage=Stage,Layer=Layer,CluterProfiler=str_trim(Genecards_index, side="left") )
	
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
boolean_counts_terms<-df_count_terms[,2:4]>0
###########################################################################################################################
# Question 1: terms in all stages
paste(names(which(boolean_counts_terms[,c("Stage_I")] & boolean_counts_terms[,c("Stage_II")] & boolean_counts_terms[,c("Stage_III")])),collapse=", ")
# Answer : Reactome:Asparagine N-linked glycosylation, Reactome:Cell surface interactions at the vascular wall, Reactome:Neutrophil degranulation, Reactome:SARS-CoV Infections
# Genes in these terms
paste(unique(df_all_annotation_per_stage[df_all_annotation_per_stage$CluterProfiler %in% names(which(boolean_counts_terms[,c("Stage_I")] & boolean_counts_terms[,c("Stage_II")] & boolean_counts_terms[,c("Stage_III")])),"Symbol"]), collapse=", ")
# "LTF, OLFM4, NPL, GRB7, TKFC, CXADR, COPB2, MOGS, RPN2, ANO1, ADAM8, ANO9, SRC, HSPA1A, ARF5, CAPN1, DERA, FUT8, AP2S1, MTA3, RPS5, SEH1L, TMED2, PSMC6, PSMD7, RPS16, AKT2, PTGES3, QPCT, EPCAM, SRD5A3, IDH1, SLC16A3, CCT8, CDA, B4GALT3, PSMC2, BRMS1"
######################################################################################################################
# Question 2: terms in Stage I only
paste(names(which(boolean_counts_terms[,c("Stage_I")] & !boolean_counts_terms[,c("Stage_II")] & !boolean_counts_terms[,c("Stage_III")])),collapse=", ")
# Answer : GO:cytoplasmic stress granule, GO:specific granule lumen, GO:tertiary granule lumen
paste(unique(df_all_annotation_per_stage[df_all_annotation_per_stage$CluterProfiler %in% names(which(boolean_counts_terms[,c("Stage_I")] & !boolean_counts_terms[,c("Stage_II")] & !boolean_counts_terms[,c("Stage_III")])),"Symbol"]), collapse=", ")
# LTF, OLFM4, GRB7, DHX36
######################################################################################################################
# Question 3: terms in Stage II only
paste(names(which(!boolean_counts_terms[,c("Stage_I")] & boolean_counts_terms[,c("Stage_II")] & !boolean_counts_terms[,c("Stage_III")])),collapse=", ")
# Answer : there is not
######################################################################################################################
# Question 4: terms in Stages II and III
paste(names(which(!boolean_counts_terms[,c("Stage_I")] & boolean_counts_terms[,c("Stage_II")] & boolean_counts_terms[,c("Stage_III")])),collapse=", ")
# Answer : Reactome:Mitotic G1 phase and G1/S transition, Reactome:Neddylation, Reactome:Cellular response to chemical stress, Reactome:Translation, Reactome:Intracellular signaling by second messengers, Reactome:PIP3 activates AKT signaling
paste(unique(df_all_annotation_per_stage[df_all_annotation_per_stage$CluterProfiler %in% names(which(!boolean_counts_terms[,c("Stage_I")] & boolean_counts_terms[,c("Stage_II")] & boolean_counts_terms[,c("Stage_III")])),"Symbol"]), collapse=", ")
# HIF1A, MRPL19, RPN2, AARS2, RPLP1, PRKCI, SRC, MTA3, RPS5, PSMC6, PSMA7, HM13, PSMD7, RPS16, AKT2, RHEB, COX6A1, MRPL18, ORC2, COX7A2L, RNF2, CDKN2C, MRPS7, UBE2M, RPA1, IDH1, SRP9, ELOC, MRPL17, PSMD4, PSMC2, NUDT2, COPS6, E2F6, TRMT112, MRPL48, RPS6KB2, GADD45GIP1, MRPS23, MRPS16, MRPL40, AKT1S1
######################################################################################################################
# Question 5: terms in Stages III only
paste(names(which(!boolean_counts_terms[,c("Stage_I")] & !boolean_counts_terms[,c("Stage_II")] & boolean_counts_terms[,c("Stage_III")])),collapse=", ")
# Answer : Reactome:PCP/CE pathway, Reactome:PTEN Regulation, Reactome:TCF dependent signaling in response to WNT, GO:ribosome, KEGG:Amyotrophic lateral sclerosis, KEGG:Huntington disease, Reactome:HIV Infection, Reactome:Signaling by WNT, KEGG:Pathways of neurodegeneration - multiple diseases, GO:mitochondrial matrix, GO:mitochondrial protein-containing complex, KEGG:Alzheimer disease, GO:mitochondrial inner membrane
terms_stage_III<-names(which(!boolean_counts_terms[,c("Stage_I")] & !boolean_counts_terms[,c("Stage_II")] & boolean_counts_terms[,c("Stage_III")]))
genes_stage_III<-unique(df_all_annotation_per_stage[df_all_annotation_per_stage$CluterProfiler %in% names(which(!boolean_counts_terms[,c("Stage_I")] & !boolean_counts_terms[,c("Stage_II")] & boolean_counts_terms[,c("Stage_III")])),"Symbol"])
# Reactome:PCP/CE pathway, Reactome:PTEN Regulation, Reactome:TCF dependent signaling in response to WNT, GO:ribosome, KEGG:Amyotrophic lateral sclerosis, KEGG:Huntington disease, Reactome:HIV Infection, Reactome:Signaling by WNT, KEGG:Pathways of neurodegeneration - multiple diseases, GO:mitochondrial matrix, GO:mitochondrial protein-containing complex, KEGG:Alzheimer disease, GO:mitochondrial inner membrane
# DVL2, PSMB1, VTA1, CAPN1, AP2S1, MTA3, RPS5, SEH1L, PSMC6, TIMM9, PSMA7, PSMD7, RPS16, AKT2, RHEB, NDUFS8, COX6A1, MRPL18, COX7A2L, BCL9, PLA2G4A, RNF2, CLTA, MRPS7, NXT1, ELOC, SUPV3L1, MALSU1, MRPL17, PSMD4, PSMC2, MAIP1, ABCE1, NUDT2, PHB1, MRPL48, GADD45GIP1, MRPS23, MRPS16, GLRX5, NELFA, MRPL40, HSPA14, SLC39A10, TIMM23
paste(terms_stage_III, collapse=", ")
paste(genes_stage_III, collapse=", ")
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

# Vertice colours of genes
#V(graph_all_stages)$labels                                                                        <- ""                                                                                                                                                                                   <-"circle"


terms_stage_III<-names(which(!boolean_counts_terms[,c("Stage_I")] & !boolean_counts_terms[,c("Stage_II")] & boolean_counts_terms[,c("Stage_III")]))
genes_stage_III<-unique(df_all_annotation_per_stage[df_all_annotation_per_stage$CluterProfiler %in% names(which(!boolean_counts_terms[,c("Stage_I")] & !boolean_counts_terms[,c("Stage_II")] & boolean_counts_terms[,c("Stage_III")])),"Symbol"])

V(graph_all_stages)[which(names(V(graph_all_stages)) %in% terms_stage_III)]$labels    <- terms_stage_III
V(graph_all_stages)[which(names(V(graph_all_stages)) %in% genes_stage_III)]$labels     <- genes_stage_III

# FindClusters_resolution
png(filename=paste(output_folder,"Plot_Stage_all_per_Stagge.png",sep=""), width = 30, height = 30, res=600, units = "cm")
	#plot(graph_all_stages, layout=layout_nicely, vertex.label=NA) # Stage I
	#plot(graph_all_stages, layout=  layout_with_kk, vertex.label=NA) # Stage II
	#plot(graph_all_stages, layout=   layout_nicely, vertex.label=V(graph_all_stages)$labels) # Stage II
	plot(graph_all_stages, layout=   layout_nicely,vertex.label=NA ) # Stage II
dev.off()

plot(graph_all_stages, layout=   layout_nicely, vertex.label=NA) # Stage II
tkplot(graph_all_stages, layout=   layout_nicely, vertex.label=NA)

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

