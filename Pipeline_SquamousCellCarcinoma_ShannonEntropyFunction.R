########################################################################################################################################
# If at least one of the genes in the pair are in the interactome
interactome_data_stage_I<-rbind(interactome_data[interactome_data$Gene1 %in% genes_interactome_stage_I,],
interactome_data[interactome_data$Gene2 %in% genes_interactome_stage_I,])

# If at least one of the genes in the pair are in the interactome
interactome_data_stage_II<-rbind(interactome_data[interactome_data$Gene1 %in% genes_interactome_stage_II,],
interactome_data[interactome_data$Gene2 %in% genes_interactome_stage_II,])

# If at least one of the genes in the pair are in the interactome
interactome_data_stage_III<-rbind(interactome_data[interactome_data$Gene1 %in% genes_interactome_stage_III,],
interactome_data[interactome_data$Gene2 %in% genes_interactome_stage_III,])
########################################################################################################################################
# copy interactome_data_stages
interactome_data_stage_I_clean<-interactome_data_stage_I
interactome_data_stage_II_clean<-interactome_data_stage_II
interactome_data_stage_III_clean<-interactome_data_stage_III

merge_interactome_data<-rbind(data.frame(Gene1=interactome_data_stage_I_clean$Gene1,Gene2=interactome_data_stage_I_clean$Gene2,Stage="Stage I"),
data.frame(Gene1=interactome_data_stage_II_clean$Gene1,Gene2=interactome_data_stage_II_clean$Gene2,Stage="Stage II"),
data.frame(Gene1=interactome_data_stage_III_clean$Gene1,Gene2=interactome_data_stage_III_clean$Gene2,Stage="Stage III"))

# Clean the tables
for (gene_pair_index in rownames(merge_interactome_data))
{
    # interactome_data_stage
    pair_gene_id_I <-merge_interactome_data[gene_pair_index,"Gene1"]
    pair_gene_id_II<-merge_interactome_data[gene_pair_index,"Gene2"]

    # Re-order gene ids
    if(pair_gene_id_II<pair_gene_id_I)
    {
      merge_interactome_data[gene_pair_index,"Gene1"]<-pair_gene_id_II
      merge_interactome_data[gene_pair_index,"Gene2"]<-pair_gene_id_I     
    }
    # Re-order gene ids
    if(pair_gene_id_II==pair_gene_id_I)
    {
      merge_interactome_data[gene_pair_index,"Gene2"]<-"REPEAT"
    }  
}
  
merge_interactome_data<-unique(merge_interactome_data)
########################################################################################################################################

########################################################################################################################################
df_stageI_connectivity   <-unique(data.frame(Conectivity=table(c(interactome_data_stage_I$Gene1,interactome_data_stage_I$Gene2))))
df_stageII_connectivity  <-unique(data.frame(Conectivity=table(c(interactome_data_stage_II$Gene1,interactome_data_stage_II$Gene2))))
df_stageIII_connectivity <-unique(data.frame(Conectivity=table(c(interactome_data_stage_III$Gene1,interactome_data_stage_III$Gene2))))
########################################################################################################################################
