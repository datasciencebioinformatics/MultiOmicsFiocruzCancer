#######################################################################################################################################
# Path to files of selected_genes                                                                                                     # 
selected_genes_Stage_I_file       <-"/home/felipe/Documentos/LungPortal/samples/pos_unique_genes_Stage_I.tsv"                         #
selected_genes_Stage_II_file      <-"/home/felipe/Documentos/LungPortal/samples/pos_unique_genes_Stage_II.tsv"                        #
selected_genes_Stage_III_file     <-"/home/felipe/Documentos/LungPortal/samples/pos_unique_genes_Stage_III.tsv"                       #
#######################################################################################################################################
# Load data                                                                                                                           #
selected_genes_Stage_I_data       <-read.table(file = selected_genes_Stage_I_file, sep = '\t', header = TRUE,fill=TRUE)               #
selected_genes_Stage_II_data      <-read.table(file = selected_genes_Stage_II_file, sep = '\t', header = TRUE,fill=TRUE)              #
selected_genes_Stage_III_data     <-read.table(file = selected_genes_Stage_III_file, sep = '\t', header = TRUE,fill=TRUE)             #
#######################################################################################################################################
library("rrvgo")
stage_I_sim_matrix <- calculateSimMatrix(selected_genes_Stage_I_data$Gene,orgdb="org.Hs.eg.db",ont="BP", method="Rel")
stage_II_sim_matrix <- calculateSimMatrix(selected_genes_Stage_II_data$Gene,orgdb="org.Hs.eg.db",ont="BP", method="Rel")
stage_III_sim_matrix <- calculateSimMatrix(selected_genes_Stage_III_data$Gene,orgdb="org.Hs.eg.db",ont="BP", method="Rel")
#######################################################################################################################################
library("topGO")
sampleGOdata <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = list(selected_genes_Stage_I_data$Gene, selected_genes_Stage_II_data$Gene, selected_genes_Stage_III_data$Gene), geneSel = list(selected_genes_Stage_I_data$Gene, selected_genes_Stage_II_data$Gene, selected_genes_Stage_III_data$Gene), nodeSize = 10,annot = annFUN.db, affyLib = affyLib)
#######################################################################################################################################
# https://bioconductor.org/packages/release/bioc/vignettes/rrvgo/inst/doc/rrvgo.html
scores <- setNames(-log10(stage_I_sim_matrix$qvalue), stage_I_sim_matrix$ID)
reducedTerms <- reduceSimMatrix(stage_I_sim_matrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
#######################################################################################################################################
