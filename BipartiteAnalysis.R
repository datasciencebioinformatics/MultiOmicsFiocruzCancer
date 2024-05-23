############################################################################################################################################
library(RCy3)
library(RColorBrewer)
# cyrestGET(operation = NULL, parameters = NULL, base.url = "http://localhost:1234")
############################################################################################################################################
coexpression_network_Stage_I   <-"/home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/df_stageI_interactome.tsv"
coexpression_network_Stage_II  <-"/home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/df_stageII_interactome.tsv"
coexpression_network_Stage_III <-"/home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/df_stageIII_interactome.tsv"
##################################################################################################
interactome_network_Stage_I   <-"/home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/df_stageI_interactome_interactome.tsv"
interactome_network_Stage_II  <-"/home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/df_stageII_interactome_interactome.tsv"
interactome_network_Stage_III <-"/home/felipe/Documentos/LungPortal/output/threhold_RPKM_3_threhold_log2foldchange_1.0_FDR_0.05_threhold_correlation_0.99/df_stageIII_interactome_interactome.tsv"
##################################################################################################
coexpression_network_Stage_I     <-read.table(file = coexpression_network_Stage_I, sep = '\t', header = TRUE,fill=TRUE)    #
coexpression_network_Stage_II     <-read.table(file = coexpression_network_Stage_II, sep = '\t', header = TRUE,fill=TRUE)    #
coexpression_network_Stage_III     <-read.table(file = coexpression_network_Stage_III, sep = '\t', header = TRUE,fill=TRUE)    #
##################################################################################################
interactome_network_Stage_I     <-read.table(file = interactome_network_Stage_I, sep = '\t', header = TRUE,fill=TRUE)    #
interactome_network_Stage_II     <-read.table(file = interactome_network_Stage_II, sep = '\t', header = TRUE,fill=TRUE)    #
interactome_network_Stage_III     <-read.table(file = interactome_network_Stage_III, sep = '\t', header = TRUE,fill=TRUE)    #
#######################################################################################################################################
# Merge data.frame
merge_interactome_gene_symbol <- merge_interactome_gene_symbol[match(unique(merge_interactome_gene_symbol$gene_id), merge_interactome_gene_symbol$gene_id),]

# Assert rownames
rownames(merge_interactome_gene_symbol)<-merge_interactome_gene_symbol$gene_id

# Converted gene_ids
coexpression_network_Stage_I    <-na.omit(data.frame(Gene1=merge_interactome_gene_symbol[coexpression_network_Stage_I$Gene1,"gene_symbol"],Gene2=merge_interactome_gene_symbol[coexpression_network_Stage_I$Gene2,"gene_symbol"]))
coexpression_network_Stage_II    <-na.omit(data.frame(Gene1=merge_interactome_gene_symbol[coexpression_network_Stage_II$Gene1,"gene_symbol"],Gene2=merge_interactome_gene_symbol[coexpression_network_Stage_II$Gene2,"gene_symbol"]))
coexpression_network_Stage_III    <-na.omit(data.frame(Gene1=merge_interactome_gene_symbol[coexpression_network_Stage_III$Gene1,"gene_symbol"],Gene2=merge_interactome_gene_symbol[coexpression_network_Stage_III$Gene2,"gene_symbol"]))

# Read igraph
interactome_igraph_stage_I<-graph_from_data_frame(stage_I_interaction_data, directed = FALSE, vertices = NULL)
interactome_igraph_stage_II<-graph_from_data_frame(stage_II_interaction_data, directed = FALSE, vertices = NULL)
interactome_igraph_stage_III<-graph_from_data_frame(stage_III_interaction_data, directed = FALSE, vertices = NULL)

# greedy method (hiearchical, fast method)
cluster_igraph_stage_I  <- cluster_fast_greedy(interactome_igraph_stage_I)
cluster_igraph_stage_II <- cluster_fast_greedy(interactome_igraph_stage_II)
cluster_igraph_stage_III <- cluster_fast_greedy(interactome_igraph_stage_III)
############################################################################################################################################
# First analysis, co-expression network and sub-interactome network for each stage, together with bipartite networks.
# Co-expression network for stage I
# Co-expression network for stage II
# Co-expression network for stage III

# Set colour
nb.cols  <- max(c(length(unique(membership(cluster_igraph_stage_I))),length(unique(membership(cluster_igraph_stage_II))), length(unique(membership(cluster_igraph_stage_III)))))
mycolors <- data.frame(colour=colorRampPalette(brewer.pal(8, "Set3"))(nb.cols))

# Generate colors based on media type:
V(interactome_igraph_stage_I)$color  <- mycolors[as.numeric( membership(cluster_igraph_stage_I)),]
V(interactome_igraph_stage_II)$color <- mycolors[as.numeric( membership(cluster_igraph_stage_II)),]
V(interactome_igraph_stage_III)$color<- mycolors[as.numeric( membership(cluster_igraph_stage_III)),]

# Set node size based on audience size:
V(interactome_igraph_stage_I)$size <- degree(interactome_igraph_stage_I)
V(interactome_igraph_stage_II)$size <- degree(interactome_igraph_stage_II)
V(interactome_igraph_stage_III)$size <- degree(interactome_igraph_stage_III)

# Sub-interactome network for stage I
# Sub-interactome network for stage II
# Sub-interactome network for stage III

# Bi-partite Co-expression vs. Sub-interactome network for stage I
# Bi-partite Co-expression vs. Sub-interactome network for stage II
# Bi-partite Co-expression vs. Sub-interactome network for stage III
################################################################################################
# Planejamento igraph networks
# 1- Primeira análise, rede de coexpressão e rede de subinteractomas para cada estágio, juntamente com redes bipartidas.
# Rede co-expressão Estágio I <-> Bipartite co-expressão com sub-interactoma Estágio I <-> Rede sub-interactoma Estágio I
# Rede co-expressão Estágio II <-> Bipartite co-expressão com sub-interactoma Estágio II <-> Rede sub-interactoma Estágio II
# Rede co-expressão Estágio III <-> Bipartite co-expressão com sub-interactoma Estágio II <-> Rede sub-interactoma Estágio III

# 2- Segunda análise, rede para cada estágio, juntamente com redes bipartidas para cada par de estágios
Rede estágio I- Bipartite estagios I e II - Rede estágio II - bipartites II e III - Rede estágio III
