library(DescTools)

# Take stageI_list_of_genes
list_of_genes<-unique(merge_interactome_gene_symbol[,c("gene_id","PPI")])

# File path
file_genes_Stage_I   <- "/home/felipe/Documentos/LungPortal/samples/log2foldchange_stageI_pos.tsv"
file_genes_Stage_II  <- "/home/felipe/Documentos/LungPortal/samples/log2foldchange_stageII_pos.tsv"
file_genes_Stage_III <- "/home/felipe/Documentos/LungPortal/samples/log2foldchange_stageIII_pos.tsv"

# Gene table
genes_Stage_I       <-read.table(file = file_genes_Stage_I, sep = '\t', header = TRUE,fill=TRUE)         #
genes_Stage_II      <-read.table(file = file_genes_Stage_II, sep = '\t', header = TRUE,fill=TRUE)#
genes_Stage_III     <-read.table(file = file_genes_Stage_III, sep = '\t', header = TRUE,fill=TRUE)                 #


genes_Stage_I<-list_of_genes[list_of_genes$gene %in% genes_Stage_I$Gene,]
genes_Stage_II<-list_of_genes[list_of_genes$gene %in% genes_Stage_II$Gene,]
genes_Stage_III<-list_of_genes[list_of_genes$gene %in% genes_Stage_III$Gene,]

round(Entropy(genes_Stage_I$PPI, base = 2),3)
round(Entropy(genes_Stage_II$PPI, base = 2),3)
round(Entropy(genes_Stage_III$PPI, base = 2),3)
