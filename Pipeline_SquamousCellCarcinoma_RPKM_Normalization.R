unstranded_data_file                <- "/home/felipe/Documentos/LungPortal/samples/unstranded.tsv"                        #
colData_file                        <- "/home/felipe/Documentos/LungPortal/samples/colData.tsv"                           #
merge_interactome_file      	    <-"/home/felipe/Documentos/LungPortal/samples/merge_interactome_gene_symbol"
unstranded_data                     <-read.table(file = unstranded_data_file, sep = '\t', header = TRUE,fill=TRUE)         #
colData_data                        <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)       
merge_interactome_gene_symbol	    <-read.table(file = merge_interactome_file, sep = '\t', header = TRUE,fill=TRUE)     #


####################################################################################################################
# RPKM normalization
# Run DESeq2
dds_stages <- DESeqDataSetFromMatrix(countData = unstranded_data, colData=colData_data, design = as.formula(paste("~ 0+ stages"))  )

# Run DESeq2
dds_stages <- DESeq(dds_stages)

print(resultsNames(dds_stages)[4])		
#####################################################################################################################
# Run varianceStabilizingTransformation                                                                             #
vst_stages_sub<-varianceStabilizingTransformation(dds_stages, blind = TRUE, fitType = "parametric")                 #
                                                                                                                    #
                                                                                                                    #
# Obtain normalized coutns                                                                                          #
norm_counts<-counts(dds_stages, normalized = TRUE)                                                      #
#####################################################################################################################
df_norm_conts<-melt(norm_counts)
###################################################################################################################
unstranded_data<-na.omit(unstranded_data[names(ebt),])

# RPKM normalization
# Run DESeq2
# Run DESeq2
dds_stages <- DESeqDataSetFromMatrix(countData = unstranded_data, colData=colData_data, design = as.formula(paste("~ 0+ stages"))  )

# Run DESeq2
dds_stages <- DESeq(dds_stages)

# Df s6tages I
df_stages<-data.frame(results(dds_stages))
#####################################################################################################################
# Run varianceStabilizingTransformation                                                                             #
vst_stages_sub<-varianceStabilizingTransformation(dds_stages, blind = TRUE, fitType = "parametric")                 #
                                                                                                                    #
                                                                                                                    #
# Obtain normalized coutns                                                                                          #
norm_counts<-counts(dds_stages, normalized = TRUE)                                                      #
#####################################################################################################################
library(AnnotationHub)
ah <- AnnotationHub()
edb <- query(ah, c("EnsDb", "Hsapiens", "v87"))[[1]]
txdf <- select(edb,keys=keys(edb, "GENEID"),columns=c("GENEID","TXID"),keytype="GENEID")
ebt <- exonsBy(edb, by="gene")

# merge_interactome_gene_sub
merge_interactome_gene_sub<-merge_interactome_gene_symbol[which(merge_interactome_gene_symbol$Gene_id %in% names(ebt)),]
merge_interactome_gene_sub<-merge_interactome_gene_symbol[which(merge_interactome_gene_symbol$gene_id %in% rownames(norm_counts)),]
df_unique_genes<-unique(data.frame(gene_id=merge_interactome_gene_sub$gene_id,Gene_id=merge_interactome_gene_sub$Gene_id))


# EBT
ebt<-ebt[names(ebt)[which(names(ebt) %in% df_unique_genes$Gene_id)]]

df_unique_genes<-unique(df_unique_genes[which(df_unique_genes$Gene_id %in% names(ebt)),c("gene_id","Gene_id")])
rownames(df_unique_genes)<-df_unique_genes$gene_id
names(ebt)<-df_unique_genes[names(ebt),"gene_id"]

# RPKM normalization
# Run DESeq2
# Run DESeq2
dds_stages <- DESeqDataSetFromMatrix(countData = na.omit(unstranded_data[names(ebt),]), colData=colData_data, design = as.formula(paste("~ 0+ stages"))  )

# Run DESeq2
dds_stages <- DESeq(dds_stages)

ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

names(ebt)<-df_unique_genes[names(ebt),"gene_id"]

# EBT
ebt<-ebt[names(ebt)[which(names(ebt) %in% rownames(dds_stages))]]

# Df s6tages I
#df_stages<-na.omit(data.frame(results(dds_stages)))

rowRanges(dds_stages) <- GRangesList(ebt)

# fpkm
dds_stages_fpkm<-fpkm(dds_stages)


df_unique_genes


df_genes<-na.omit(df_unique_genes[rownames(unstranded_data),])
library (edgeR)
library (EDASeq)
library("biomaRt")
ensembl_list <- c("ENSG00000000003","ENSG00000000419","ENSG00000000457","ENSG00000000460")
geneLengthAndGCContent_1<-getGeneLengthAndGCContent(df_genes$Gene_id[1:1000], "hsa")

rownames(df_unique_genes)<-df_unique_genes$Gene_id
df_unique_genes[rownames(geneLengthAndGCContent_1),"gene_id"]

rpkm(unstranded_data[df_unique_genes[rownames(geneLengthAndGCContent_1),"gene_id"],], gene.length = data.frame(geneLengthAndGCContent_1)$length)

ggplot(df_stages, aes(log2FoldChange)) +
  geom_histogram(color = "#000000", fill = "#0099F8",bin=10) +
  geom_vline(aes(xintercept = mean(value)), color = "#000000", size = 1.25) +
  geom_vline(aes(xintercept = mean(value) + var(log2FoldChange)), color = "#000000", size = 1, linetype = "dashed") +      ggtitle(paste("Stage I vs. Stages II and III\nlog2FoldChange"))+ theme(legend.position='bottom')
