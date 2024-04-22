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
# RPKM normalization
# The normalization is done for each 1000 genes duo to limitation in the biomart connection

# First, gene length and gc content for all genes in the reads count table
# Take the gene names, without variant identification
# vector to store all gene ids
gene_ids<-c()
for (gene_id in rownames(unstranded_data)) 
{
    # Store gene ids
    gene_ids<-c(gene_ids,strsplit(gene_id,".",fixed=T)[[1]][[1]])
}
# Split gene_ids vector in parts
gene_ids_vector<-split(gene_ids,ceiling(seq_along(gene_ids) / 1000))

# Data.frame to store geneLengthAndGCContent
df_geneLengthAndGCContent<-data.frame(length=c(),gc=c())

# For each part of the vectors
for (index in names(gene_ids_vector) )
{    
    # Concatenate files
    df_geneLengthAndGCContent<-rbind(df_geneLengthAndGCContent,getGeneLengthAndGCContent(gene_ids_vector[[index]], "hsa"))
}
####################################################################################################################
unstranded_rpkm<-rpkm(unstranded_data[df_unique_genes[rownames(geneLengthAndGCContent_1),"gene_id"],], gene.length = data.frame(geneLengthAndGCContent_1)$length)
####################################################################################################################
# Save TSV file with genes from Stage1
write_tsv(unstranded_rpkm, paste(output_dir,"DESeq2_selected_genes_Stage_pos_",stage_index,".tsv",sep=""))
