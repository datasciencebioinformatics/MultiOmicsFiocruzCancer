#######################################################################################################################################
# Path to files of selected_genes                                                                                                     # 
selected_genes_Stage_I_file       <-"/home/felipe/Documentos/LungPortal/samples/pos_unique_genes_Stage_I.tsv"                         #
selected_genes_Stage_II_file      <-"/home/felipe/Documentos/LungPortal/samples/pos_unique_genes_Stage_II.tsv"                        #
selected_genes_Stage_III_file     <-"/home/felipe/Documentos/LungPortal/samples/pos_unique_genes_Stage_III.tsv"                       #
unstranded_file                   <- "/home/felipe/Documentos/LungPortal/samples/unstranded.rna_seq.augmented_star_gene_counts.tsv"   #
#######################################################################################################################################
# Load data                                                                                                                           #
selected_genes_Stage_I_data       <-read.table(file = selected_genes_Stage_I_file, sep = '\t', header = TRUE,fill=TRUE)               #
selected_genes_Stage_II_data      <-read.table(file = selected_genes_Stage_II_file, sep = '\t', header = TRUE,fill=TRUE)              #
selected_genes_Stage_III_data     <-read.table(file = selected_genes_Stage_III_file, sep = '\t', header = TRUE,fill=TRUE)             #
#######################################################################################################################################
# Run DESeq2
dds_stages <- DESeqDataSetFromMatrix(countData = unstranded_data, colData=colData[colnames(unstranded_data),], design = ~  age_at_index +  gender +stages)

# Run DESeq2
dds_stages <- DESeq(dds_stages)

# Obtain differential expression numbers
resultsNames(dds_stages)

# Df s6tages I
df_stages<-data.frame(results(dds_stages,name="stages_Stage.II_vs_Stage.I"))
#######################################################################################################################################
library("rrvgo")
stage_I_sim_matrix <- calculateSimMatrix(selected_genes_Stage_I_data$Gene,orgdb="org.Hs.eg.db",ont="BP", method="Rel")
stage_II_sim_matrix <- calculateSimMatrix(selected_genes_Stage_II_data$Gene,orgdb="org.Hs.eg.db",ont="BP", method="Rel")
stage_III_sim_matrix <- calculateSimMatrix(selected_genes_Stage_III_data$Gene,orgdb="org.Hs.eg.db",ont="BP", method="Rel")
#######################################################################################################################################
# Fist take 
#######################################################################################################################################
library("topGO")
library("hgu95av2.db")

# Vectors for the pvalue of selected genes
vector_Stage_I <-selected_genes_Stage_I_data$pvalue
vector_Stage_II <-selected_genes_Stage_II_data$pvalue
vector_Stage_III <-selected_genes_Stage_III_data$pvalue
vector_all <- df_stages$pvalue

# Names of genes
names(vector_Stage_I)<-selected_genes_Stage_I_data$Gene
names(vector_Stage_II)<-selected_genes_Stage_II_data$Gene
names(vector_Stage_III)<-selected_genes_Stage_III_data$Gene
names(vector_all) <- rownames(df_stages)

## To do:
## Check how implement function to consider folchange and padj threshold here
## Otherwise, check if disconsider positive and negatives
## Otherwise, check if disconsider foldchange values
# Question is : I can use a vector, but have one parameter
# Consideration is : I am using this analysis to check biological relevance of signals
# I can use all DE genes padj < 0.05

## create a gene selection function to select significant genes
topDiffGenes <- function(padj) {return (padj < 0.05)}

# Check how to use Ensemble transcripts ID in topgo
topGO = new("topGOdata", description="stages", ontology= "BP",  allGenes = vector_Stage_I, geneSel = topDiffGenes, nodeSize = 10, annot=annFUN.org, mapping="org.Hs.eg.db", ID = "GeneName")

sampleGOdata <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = geneList, geneSel = topDiffGenes, nodeSize = 10,annot = annFUN.db, affyLib = affyLib)
#######################################################################################################################################
# https://bioconductor.org/packages/release/bioc/vignettes/rrvgo/inst/doc/rrvgo.html
scores <- setNames(-log10(stage_I_sim_matrix$qvalue), stage_I_sim_matrix$ID)
reducedTerms <- reduceSimMatrix(stage_I_sim_matrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
#######################################################################################################################################
