#######################################################################################################################
hsa <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://asia.ensembl.org")#
#######################################################################################################################
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
####################################################################################################################
# Get gene coordinates
gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), filters="ensembl_gene_id", values=df_gene_ids$gene_id, mart=hsa)

# Get gene size
gene_coords$size=gene_coords$end_position - gene_coords$start_position

# Set colnames
colnames(df_gene_ids)<-c("gene_id","transcript_id")
colnames(gene_coords)[2]<-"gene_id"

# Merge data
merged_data<-merge(df_gene_ids,gene_coords,by="gene_id")

# Calculate 
unstranded_rpkm<-rpkm(unstranded_data[merged_data$transcript_id,], gene.length = merged_data$size) #
####################################################################################################################
write.table(data.frame(unstranded_rpkm), paste(output_dir,"unstranded_rpkm.tsv",sep="/"), sep = '\t')
