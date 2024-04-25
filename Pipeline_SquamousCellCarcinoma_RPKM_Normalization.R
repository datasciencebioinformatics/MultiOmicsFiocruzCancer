library(AnnotationHub)
library (edgeR)
library (EDASeq)
library("biomaRt")
####################################################################################################################
# A script to normalize reads count to RPKM                                                                        #
####################################################################################################################
# Path to file                                                                                                     #
unstranded_data_file                <- "/home/felipe/Documentos/LungPortal/samples/unstranded.tsv"                 #
colData_file                        <- "/home/felipe/Documentos/LungPortal/samples/colData.tsv"                    #
merge_interactome_file      	    <-"/home/felipe/Documentos/LungPortal/samples/merge_interactome_gene_symbol"     #
####################################################################################################################
# Load the files                                                                                                   #
unstranded_data                     <-read.table(file = unstranded_data_file, sep = '\t', header = TRUE,fill=TRUE) #
colData_data                        <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)         #
merge_interactome_gene_symbol	      <-read.table(file = merge_interactome_file, sep = '\t', header = TRUE,fill=TRUE) # 
####################################################################################################################
# RPKM normalization                                                                                               #
# The normalization is done for each 1000 genes duo to limitation in the biomart connection                        #
# First, gene length and gc content for all genes in the reads count table                                         #
# Take the gene names, without variant identification                                                              #
# vector to store all gene ids                                                                                     #
df_gene_ids<-data.frame(gene_id=c(),gene_id_cp=c())                                                                #
for (gene_id in rownames(unstranded_data))                                                                         #
{                                                                                                                  #
    # Store gene ids                                                                                               #
    print(gene_id)
    gene_ids<-strsplit(gene_id,".",fixed=T)[[1]][[1]]                                                              #                                                  
                                                                                                                   #
    # Contatenate gene lists                                                                                       #
    df_gene_ids<-rbind(df_gene_ids,data.frame(gene_id=gene_ids,gene_id_cp=gene_id))                                #
}                                                                                                                  #
# Split gene_ids vector in parts                                                                                   #
gene_ids_vector<-split(df_gene_ids$gene_id,ceiling(seq_along(df_gene_ids$gene_id) / 1000))                         #
####################################################################################################################
# Data.frame to store geneLengthAndGCContent                                                                       #
df_geneLengthAndGCContent<-data.frame(length=c(),gc=c())                                                           #
                                                                                                                   #
# For each part of the vectors                                                                                     #
for (index in names(gene_ids_vector) )                                                                             #
{                                                                                                                  #
    # Concatenate files                                                                                            ########
    df_geneLengthAndGCContent<-rbind(df_geneLengthAndGCContent,getGeneLengthAndGCContent(gene_ids_vector[[index]], "hsa"))#
}                                                                                                                         #
rownames(df_geneLengthAndGCContent)[!grepl(".", rownames(df_geneLengthAndGCContent), fixed=TRUE)]                         #
####################################################################################################################################################################
unstranded_rpkm<-rpkm(unstranded_data[df_gene_ids$gene_id_cp,], gene.length = data.frame(df_geneLengthAndGCContent)$length) #
###################################################################################################################################################################
# Filter genes by RPKM, rowMeans(unstranded_rpkm)>5 
unstranded_rpkm<-unstranded_rpkm[rowMeans(unstranded_rpkm)>3,]
###################################################################################################################################################################

# Save TSV file with genes from Stage1                                                                                                                             #
write.table(unstranded_rpkm, file = "/home/felipe/Documentos/LungPortal/samples/unstranded_rpkm.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, append=FALSE)                  #
####################################################################################################################################################################
