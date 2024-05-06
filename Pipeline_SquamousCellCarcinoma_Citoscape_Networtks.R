############################################################################################################################################
library(RCy3)
# cyrestGET(operation = NULL, parameters = NULL, base.url = "http://localhost:1234")
############################################################################################################################################
output_folder<-"/home/felipe/Documentos/LungPortal/output/modules/"
############################################################################################################################################
df_correlation_net_stage_I<-data.frame(na.omit(unstranded_data_filter[genes_Stage_I$gene,]))
df_correlation_net_stage_II<-data.frame(na.omit(unstranded_data_filter[genes_Stage_II$gene,]))
df_correlation_net_stage_III<-data.frame(na.omit(unstranded_data_filter[genes_Stage_III$gene,]))
#######################################################################################################################################
rownames(interactome_data_stage_I)
rownames(interactome_data_stage_I)
rownames(interactome_data_stage_I)

# Merge data.frame
merge_interactome_gene_symbol <- merge_interactome_gene_symbol[match(unique(merge_interactome_gene_symbol$gene_id), merge_interactome_gene_symbol$gene_id),]

# Assert rownames
rownames(merge_interactome_gene_symbol)<-merge_interactome_gene_symbol$gene_id

# Converted gene_ids
interactome_smbols_stage_I    <-na.omit(data.frame(Gene1=merge_interactome_gene_symbol[interactome_data_stage_I$Gene1,"gene_symbol"],Gene2=merge_interactome_gene_symbol[interactome_data_stage_I$Gene2,"gene_symbol"]))
interactome_smbols_stage_II   <-na.omit(data.frame(Gene1=merge_interactome_gene_symbol[interactome_data_stage_II$Gene1,"gene_symbol"],Gene2=merge_interactome_gene_symbol[interactome_data_stage_II$Gene2,"gene_symbol"]))
interactome_smbols_stage_III  <-na.omit(data.frame(Gene1=merge_interactome_gene_symbol[interactome_data_stage_III$Gene1,"gene_symbol"],Gene2=merge_interactome_gene_symbol[interactome_data_stage_III$Gene2,"gene_symbol"]))

# Save TSV file with genes from Stage1
write_tsv(interactome_data_stage_I, paste(output_dir,"interactome_smbols_stage_I",".tsv",sep=""))
write_tsv(interactome_data_stage_II, paste(output_dir,"interactome_smbols_stage_II",".tsv",sep=""))
write_tsv(interactome_data_stage_III, paste(output_dir,"interactome_smbols_stage_III",".tsv",sep=""))
