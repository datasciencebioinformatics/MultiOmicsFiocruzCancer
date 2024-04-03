#######################################################################################################################################
# To DO:
# Here 1 : try extremes of padj threshold range (0.001 - 0.05)
# Here 2 : try extremes of topNodes 10-100
# Here 3 : check if KS classiKS can be used as score instead of fisher qvalue
# Here 4 : change threshold parameter
# Important is to improve legibility of plot
#######################################################################################################################################
# Path to files of selected_genes                                                                                                     # 
selected_genes_Stage_I_file       <-"/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_pos_stage_I.tsv"                  #
selected_genes_Stage_II_file      <-"/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_pos_stage_II.tsv"                 #
selected_genes_Stage_III_file     <-"/home/felipe/Documentos/LungPortal/output/selected_genes_Stage_pos_stage_III.tsv"                #
unstranded_file                   <-"/home/felipe/Documentos/LungPortal/samples/unstranded.rna_seq.augmented_star_gene_counts.tsv"    #
#######################################################################################################################################
# Load data                                                                                                                           #
selected_genes_Stage_I_data       <-read.table(file = selected_genes_Stage_I_file, sep = '\t', header = TRUE,fill=TRUE)               #
selected_genes_Stage_II_data      <-read.table(file = selected_genes_Stage_II_file, sep = '\t', header = TRUE,fill=TRUE)              #
selected_genes_Stage_III_data     <-read.table(file = selected_genes_Stage_III_file, sep = '\t', header = TRUE,fill=TRUE)             #
unstranded_data                   <-read.table(file = selected_genes_Stage_III_file, sep = '\t', header = TRUE,fill=TRUE)             #
#######################################################################################################################################
# Run DESeq2
dds_stages <- DESeqDataSetFromMatrix(countData = unstranded_data, colData=colData[colnames(unstranded_data),], design = ~  age_at_index +  gender +stages)

# Run DESeq2
dds_stages <- DESeq(dds_stages)

# Obtain differential expression numbers
resultsNames(dds_stages)

# Df s6tages I
df_stages<-data.frame(results(dds_stages,name=resultsNames(dds_stages)[4]))
#######################################################################################################################################
library("rrvgo")
stage_I_sim_matrix <- calculateSimMatrix(selected_genes_Stage_I_data$Gene,orgdb="org.Hs.eg.db",ont="BP", method="Rel", keytype = "ENSEMBL")
stage_II_sim_matrix <- calculateSimMatrix(substr(selected_genes_Stage_II_data$Gene, 1, 15),orgdb="org.Hs.eg.db",ont="BP", method="Rel")
stage_III_sim_matrix <- calculateSimMatrix(ssubstr(selected_genes_Stage_III_data$Gene, 1, 15),orgdb="org.Hs.eg.db",ont="BP", method="Rel")
#######################################################################################################################################
# Second, use rrvgo to plot the go terms and reduced go tems
# Try dotplot ggplot2 for comparative Stages I, II, III
# https://forum.posit.co/t/trying-to-plot-go-in-ggplot/126952
#######################################################################################################################################

library("topGO")
library("hgu95av2.db")

# Vectors for the pvalue of selected genes
vector_Stage_I <-selected_genes_Stage_I_data$padj
vector_Stage_II <-selected_genes_Stage_II_data$padj
vector_Stage_III <-selected_genes_Stage_III_data$padj
vector_all <- df_stages$pvalue

# Names of genes
names(vector_Stage_I)<-gsub("\\.\\d+", "", selected_genes_Stage_I_data$Gene)
names(vector_Stage_II)<-gsub("\\.\\d+", "", selected_genes_Stage_II_data$Gene)
names(vector_Stage_III)<-gsub("\\.\\d+", "", selected_genes_Stage_III_data$Gene)
names(vector_all) <- rownames(df_stages)

## To do:
## Check how implement function to consider folchange and padj threshold here
## Otherwise, check if disconsider positive and negatives
## Otherwise, check if disconsider foldchange values
# Question is : I can use a vector, but have one parameter
# Consideration is : I am using this analysis to check biological relevance of signals
# I can use all DE genes padj < 0.05

## create a gene selection function to select significant genes
# Here 1 : try extremes of padj threshold range (0.001 - 0.05)
topDiffGenes <- function(padj) {return (padj < 0.05)}

# Check how to use Ensemble transcripts ID in topgo
topGO_vector_Stage_I   = new("topGOdata", description="stages", ontology= "BP",  allGenes = vector_Stage_I,   geneSel = topDiffGenes, nodeSize = 10, annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ENSEMBL")
topGO_vector_Stage_II  = new("topGOdata", description="stages", ontology= "BP",  allGenes = vector_Stage_II,  geneSel = topDiffGenes, nodeSize = 10, annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ENSEMBL")
topGO_vector_Stage_III = new("topGOdata", description="stages", ontology= "BP",  allGenes = vector_Stage_III, geneSel = topDiffGenes, nodeSize = 10, annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ENSEMBL")
#######################################################################################################################################
result_topGO_vector_Stage_I <- runTest(topGO_vector_Stage_I, algorithm = "classic", statistic = "ks") # statistic = "fisher"
result_topGO_vector_Stage_II <- runTest(topGO_vector_Stage_II, algorithm = "classic", statistic = "ks") # statistic = "fisher"
result_topGO_vector_Stage_III <- runTest(topGO_vector_Stage_III, algorithm = "classic", statistic = "ks") # statistic = "fisher"

# Here 2 : try extremes of topNodes 10-100
table_topGO_vector_Stage_I   <- GenTable(topGO_vector_Stage_I,   classicKS = result_topGO_vector_Stage_I, topNodes = 10)
table_topGO_vector_Stage_II  <- GenTable(topGO_vector_Stage_II,  classicKS = result_topGO_vector_Stage_II, topNodes = 25)
table_topGO_vector_Stage_III <- GenTable(topGO_vector_Stage_III, classicKS = result_topGO_vector_Stage_III, topNodes = 25)

# https://bioconductor.org/packages/release/bioc/vignettes/rrvgo/inst/doc/rrvgo.html
simMatrix_stage_I <- calculateSimMatrix(table_topGO_vector_Stage_I$GO.ID, orgdb="org.Hs.eg.db",ont="BP",method="Rel")
simMatrix_stage_II <- calculateSimMatrix(table_topGO_vector_Stage_II$GO.ID, orgdb="org.Hs.eg.db",ont="BP",method="Rel")
simMatrix_stage_III <- calculateSimMatrix(table_topGO_vector_Stage_III$GO.ID, orgdb="org.Hs.eg.db",ont="BP",method="Rel")

# Here 3 : check if KS classiKS can be used as score instead of fisher qvalue
# Scores of Stage I
scores_stage_I    <- setNames(-log10(as.numeric(table_topGO_vector_Stage_I$classicKS)), table_topGO_vector_Stage_I$GO.ID)
scores_stage_II   <- setNames(-log10(as.numeric(table_topGO_vector_Stage_II$classicKS)), table_topGO_vector_Stage_II$GO.ID)
scores_stage_III  <- setNames(-log10(as.numeric(table_topGO_vector_Stage_III$classicKS)), table_topGO_vector_Stage_III$GO.ID)

# Here 4 : change threshold parameter
# Compute reduce terms
reducedTerms_stage_I <- reduceSimMatrix(simMatrix_stage_I, scores_stage_I, threshold=0.7,  orgdb="org.Hs.eg.db")
reducedTerms_stage_II <- reduceSimMatrix(simMatrix_stage_II, scores_stage_II, threshold=0.7,  orgdb="org.Hs.eg.db")
reducedTerms_stage_III <- reduceSimMatrix(simMatrix_stage_III, scores_stage_III, threshold=0.7,  orgdb="org.Hs.eg.db")
#######################################################################################################################################
treemapPlot(reducedTerms_stage_I)
treemapPlot(reducedTerms_stage_II)
treemapPlot(reducedTerms_stage_III)
#######################################################################################################################################
# FindClusters_resolution
png(filename=paste(output_dir,"reducedTerms_stage_I.png",sep=""), width = 16, height = 14, res=600, units = "cm")
	treemapPlot(reducedTerms_stage_I)
dev.off()
# FindClusters_resolution
png(filename=paste(output_dir,"reducedTerms_stage_II.png",sep=""), width = 16, height = 14, res=600, units = "cm")
	treemapPlot(reducedTerms_stage_II)
dev.off()
# FindClusters_resolution
png(filename=paste(output_dir,"reducedTerms_stage_III.png",sep=""), width = 16, height = 14, res=600, units = "cm")
	treemapPlot(reducedTerms_stage_III)
dev.off()
#######################################################################################################################################
reducedTerms_stage_I$Stage<-"Stage I"
reducedTerms_stage_II$Stage<-"Stage II"
reducedTerms_stage_III$Stage<-"Stage III"

reducedTerms_stages<-rbind(reducedTerms_stage_I,reducedTerms_stage_II)
reducedTerms_stages<-rbind(reducedTerms_stages,reducedTerms_stage_III)
#######################################################################################################################################
# plot: dot plot
# FindClusters_resolution
png(filename=paste(output_dir,"GOTerms_reducedTerms.png",sep=""), width = 20, height = 20, res=600, units = "cm")
	ggplot(data = reducedTerms_stages, aes(y = parentTerm, x = Stage, color = score)) +  geom_point(size=6) + scale_color_gradient(low = "red", high = "blue") +  theme_bw() +   ylab("") +  xlab("") +   ggtitle("GO enrichment analysis for Stages I, II and III") + scale_alpha(guide = "none")
dev.off()
#######################################################################################################################################
