library(pathview)

mypathway_Apoptosis<-"08403"
mypathway_WNT<-"04310"

genes<-c(genes_unique_Stage_I$gene_id,genes_unique_Stage_II$gene_id,genes_unique_Stage_III$gene_id)

# For each gene, add gene_id
for (gene_row in rownames(log2change_tumor_control))
{	
	# Store gene id in the vector
	# Simply trim the gene id before the "." to save it in the ENSEML format
	log2change_tumor_control[gene_row,"ENSEMBL"]<-strsplit(log2change_tumor_control[gene_row,"gene"], split = "\\.")[[1]][1]	
}

logFC<-log2change_tumor_control[log2change_tumor_control$ENSEMBL %in% as.vector(genes),"log2change"]
names(logFC)<-

ids_stage_I      <-bitr(log2change_tumor_control$ENSEMBL, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")

names(logFC)<-ids_stage_I[ids_stage_I$ENSEMBL %in% genes,"ENTREZID"]

pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_Apoptosis)
pathview(gene.data=logFC*100,species="hsa",pathway=mypathway_WNT)

