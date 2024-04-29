
library(corrr)
library(tidyverse)

# Take normal and tumor samples
normal_samples<-paired_sample_df$normal
tumor_samples<-paired_sample_df$tumor

# Take genes from each stage 
genes_stage_I<-na.omit(genes_Stage_I[genes_Stage_I$gene %in% rownames(unstranded_data_filter),"gene"])
genes_stage_II<-na.omit(genes_Stage_II[genes_Stage_II$gene %in% rownames(unstranded_data_filter),"gene"])
genes_stage_III<-na.omit(genes_Stage_III[genes_Stage_III$gene %in% rownames(unstranded_data_filter),"gene"])


# For each case, I take the sample for tumor and control
for (case in paired_sample_df$case)
{
    normal_sample<-paired_sample_df[,"normal"]
    tumor_sample<-paired_sample_df[,"tumor"]
    
    # Take genes from normal
    eprx_stage_I_normal  <-na.omit(unstranded_data_filter[genes_stage_I,normal_samples])
    eprx_stage_II_normal <-na.omit(unstranded_data_filter[genes_stage_II,normal_samples])
    eprx_stage_III_normal<-na.omit(unstranded_data_filter[genes_stage_III,normal_samples])
    
    # Take genes from normal
    eprx_stage_I_tumor  <-na.omit(unstranded_data_filter[genes_stage_I,tumor_samples])
    eprx_stage_II_tumor <-na.omit(unstranded_data_filter[genes_stage_II,tumor_samples])
    eprx_stage_III_tumor<-na.omit(unstranded_data_filter[genes_stage_III,tumor_samples])  
    
    # Calcualte correlation network normal
    eprx_stage_I_normal.cor <- data.frame(eprx_stage_I_normal %>%     t() %>% correlate() %>%     shave(upper = TRUE) %>%     stretch(na.rm = TRUE) %>%    filter(r >= 0.6))        
    eprx_stage_II_normal.cor <- data.frame(eprx_stage_II_normal %>%     t() %>% correlate() %>%    shave(upper = TRUE) %>%    stretch(na.rm = TRUE) %>%       filter(r >= 0.6))        
    eprx_stage_III_normal.cor <- data.frame(eprx_stage_III_normal %>%    t() %>% correlate() %>%    shave(upper = TRUE) %>%    stretch(na.rm = TRUE) %>%       filter(r >= 0.6))        
    
    # Calcualte correlation network tumor
    eprx_stage_I_tumor.cor <- data.frame(eprx_stage_I_tumor %>%    t() %>% correlate() %>%    shave(upper = TRUE) %>%    stretch(na.rm = TRUE) %>%       filter(r >= 0.6))        
    eprx_stage_II_tumor.cor <- data.frame(eprx_stage_II_tumor %>%   t() %>% correlate() %>%    shave(upper = TRUE) %>%   stretch(na.rm = TRUE) %>%         filter(r >= 0.6))        
    eprx_stage_III_tumor.cor <- data.frame(eprx_stage_III_tumor %>%   t() %>% correlate() %>%   shave(upper = TRUE) %>%    stretch(na.rm = TRUE) %>%        filter(r >= 0.6))        
    
    # Rename collumns normal
    colnames(eprx_stage_I_normal.cor)[1:2]<-c("Gene1","Gene2")
    colnames(eprx_stage_II_normal.cor)[1:2]<-c("Gene1","Gene2")
    colnames(eprx_stage_III_normal.cor)[1:2]<-c("Gene1","Gene2")
    
    # Rename collumns tumor
    colnames(eprx_stage_I_tumor.cor)[1:2]<-c("Gene1","Gene2")
    colnames(eprx_stage_II_tumor.cor)[1:2]<-c("Gene1","Gene2")
    colnames(eprx_stage_III_tumor.cor)[1:2]<-c("Gene1","Gene2")
    
    # calculate counts normal
    normal_stageI_connectivity   <-unique(data.frame(Conectivity=table(c(eprx_stage_I_normal.cor$Gene1,eprx_stage_I_normal.cor$Gene2))))
    normal_stageII_connectivity   <-unique(data.frame(Conectivity=table(c(eprx_stage_II_normal.cor$Gene1,eprx_stage_II_normal.cor$Gene2))))
    normal_stageIII_connectivity  <-unique(data.frame(Conectivity=table(c(eprx_stage_III_normal.cor$Gene1,eprx_stage_III_normal.cor$Gene2))))
    
    # calculate counts normal
    tumor_stageI_connectivity   <-unique(data.frame(Conectivity=table(c(eprx_stage_I_tumor.cor$Gene1,eprx_stage_I_tumor.cor$Gene2))))
    tumor_stageII_connectivity   <-unique(data.frame(Conectivity=table(c(eprx_stage_II_tumor.cor$Gene1,eprx_stage_II_tumor.cor$Gene2))))
    tumor_stageIII_connectivity  <-unique(data.frame(Conectivity=table(c(eprx_stage_III_tumor.cor$Gene1,eprx_stage_III_tumor.cor$Gene2))))
    
    # Rename collumns
    colnames(normal_stageI_connectivity)<-c("Gene","Conectivity")
    colnames(normal_stageII_connectivity)<-c("Gene","Conectivity")
    colnames(normal_stageIII_connectivity)<-c("Gene","Conectivity")
    
    # Rename collumns
    colnames(tumor_stageI_connectivity)<-c("Gene","Conectivity")
    colnames(tumor_stageII_connectivity)<-c("Gene","Conectivity")
    colnames(tumor_stageIII_connectivity)<-c("Gene","Conectivity")
    
    # Table for the calculation of entropy
    df_entropy_stage_I_normal   <-data.frame(table(normal_stageI_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
    df_entropy_stage_II_normal   <-data.frame(table(normal_stageII_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
    df_entropy_stage_III_normal   <-data.frame(table(normal_stageIII_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
    
    # Table for the calculation of entropy
    df_entropy_stage_I_tumor   <-data.frame(table(tumor_stageI_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
    df_entropy_stage_II_tumor   <-data.frame(table(tumor_stageII_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
    df_entropy_stage_III_tumor   <-data.frame(table(tumor_stageIII_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
    
    # Rename colnames
    colnames(df_entropy_stage_I_normal)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
    colnames(df_entropy_stage_II_normal)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
    colnames(df_entropy_stage_III_normal)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
    
    # Rename colnames
    colnames(df_entropy_stage_I_tumor)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
    colnames(df_entropy_stage_II_tumor)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
    colnames(df_entropy_stage_III_tumor)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
    
    # Calculate p(k)
    df_entropy_stage_I_normal$p_k<-df_entropy_stage_I_normal$count/sum(df_entropy_stage_I_normal$count)
    df_entropy_stage_II_normal$p_k<-df_entropy_stage_II_normal$count/sum(df_entropy_stage_II_normal$count)
    df_entropy_stage_III_normal$p_k<-df_entropy_stage_III_normal$count/sum(df_entropy_stage_III_normal$count)
    
    # Calculate p(k)
    df_entropy_stage_I_tumor$p_k<-df_entropy_stage_I_tumor$count/sum(df_entropy_stage_I_tumor$count)
    df_entropy_stage_II_tumor$p_k<-df_entropy_stage_II_tumor$count/sum(df_entropy_stage_II_tumor$count)
    df_entropy_stage_III_tumor$p_k<-df_entropy_stage_III_tumor$count/sum(df_entropy_stage_III_tumor$count)
    
    # Calculate log2(p(k))
    df_entropy_stage_I_normal$log2_pk<-log(df_entropy_stage_I_normal$p_k,2)
    df_entropy_stage_II_normal$log2_pk<-log(df_entropy_stage_II_normal$p_k,2)
    df_entropy_stage_III_normal$log2_pk<-log(df_entropy_stage_III_normal$p_k,2)
    
    df_entropy_stage_I_tumor$log2_pk<-log(df_entropy_stage_I_tumor$p_k,2)
    df_entropy_stage_II_tumor$log2_pk<-log(df_entropy_stage_II_tumor$p_k,2)
    df_entropy_stage_III_tumor$log2_pk<-log(df_entropy_stage_III_tumor$p_k,2)
    
    # Calculate p(k)*log2(p(k))
    df_entropy_stage_I_normal$p_k_mult_log2_pk<-df_entropy_stage_I_normal$p_k*df_entropy_stage_I_normal$log2_pk
    df_entropy_stage_II_normal$p_k_mult_log2_pk<-df_entropy_stage_II_normal$p_k*df_entropy_stage_II_normal$log2_pk
    df_entropy_stage_III_normal$p_k_mult_log2_pk<-df_entropy_stage_III_normal$p_k*df_entropy_stage_III_normal$log2_pk
    
    df_entropy_stage_I_tumor$p_k_mult_log2_pk<-df_entropy_stage_I_tumor$p_k*df_entropy_stage_I_tumor$log2_pk
    df_entropy_stage_II_tumor$p_k_mult_log2_pk<-df_entropy_stage_II_tumor$p_k*df_entropy_stage_II_tumor$log2_pk
    df_entropy_stage_III_tumor$p_k_mult_log2_pk<-df_entropy_stage_III_tumor$p_k*df_entropy_stage_III_tumor$log2_pk
    
    # Caclulate entropy value
    Entropy_stage_I_value_Carels_normal  <-abs(sum(df_entropy_stage_I_normal$p_k_mult_log2_pk))
    Entropy_stage_II_value_Carels_normal  <-abs(sum(df_entropy_stage_II_normal$p_k_mult_log2_pk))
    Entropy_stage_III_value_Carels_normal  <-abs(sum(df_entropy_stage_III_normal$p_k_mult_log2_pk))
    
    # Caclulate entropy value
    Entropy_stage_I_value_Carels_tumor  <-abs(sum(df_entropy_stage_I_tumor$p_k_mult_log2_pk))
    Entropy_stage_II_value_Carels_tumor  <-abs(sum(df_entropy_stage_II_tumor$p_k_mult_log2_pk))
    Entropy_stage_III_value_Carels_tumor  <-abs(sum(df_entropy_stage_III_tumor$p_k_mult_log2_pk))
    
    # Rename collumns
    colnames(eprx_stage_I_tumor.cor)[1:2]<-c("Gene1","Gene2")
    colnames(eprx_stage_II_tumor.cor)[1:2]<-c("Gene1","Gene2")
    colnames(eprx_stage_III_tumor.cor)[1:2]<-c("Gene1","Gene2")
    
    # calculate counts
    tumor_stageI_connectivity   <-unique(data.frame(Conectivity=table(c(eprx_stage_I_tumor.cor$Gene1,eprx_stage_I_tumor.cor$Gene2))))
    tumor_stageII_connectivity   <-unique(data.frame(Conectivity=table(c(eprx_stage_II_tumor.cor$Gene1,eprx_stage_II_tumor.cor$Gene2))))
    tumor_stageIII_connectivity  <-unique(data.frame(Conectivity=table(c(eprx_stage_III_tumor.cor$Gene1,eprx_stage_III_tumor.cor$Gene2))))
    
    # Rename collumns
    colnames(tumor_stageI_connectivity)<-c("Gene","Conectivity")
    colnames(tumor_stageII_connectivity)<-c("Gene","Conectivity")
    colnames(tumor_stageIII_connectivity)<-c("Gene","Conectivity")
    
    # Table for the calculation of entropy
    df_entropy_stage_I_tumor   <-data.frame(table(tumor_stageI_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
    df_entropy_stage_II_tumor   <-data.frame(table(tumor_stageII_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
    df_entropy_stage_III_tumor   <-data.frame(table(tumor_stageIII_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
    
    # Rename colnames
    colnames(df_entropy_stage_I_tumor)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
    colnames(df_entropy_stage_II_tumor)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
    colnames(df_entropy_stage_III_tumor)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
    
    # Calculate p(k)
    df_entropy_stage_I_tumor$p_k<-df_entropy_stage_I_tumor$count/sum(df_entropy_stage_I_tumor$count)
    df_entropy_stage_II_tumor$p_k<-df_entropy_stage_II_tumor$count/sum(df_entropy_stage_II_tumor$count)
    df_entropy_stage_III_tumor$p_k<-df_entropy_stage_III_tumor$count/sum(df_entropy_stage_III_tumor$count)
    
    # Calculate log2(p(k))
    df_entropy_stage_I_tumor$log2_pk<-log(df_entropy_stage_I_tumor$p_k,2)
    df_entropy_stage_II_tumor$log2_pk<-log(df_entropy_stage_II_tumor$p_k,2)
    df_entropy_stage_III_tumor$log2_pk<-log(df_entropy_stage_III_tumor$p_k,2)
    
    # Calculate p(k)*log2(p(k))
    df_entropy_stage_I_tumor$p_k_mult_log2_pk<-df_entropy_stage_I_tumor$p_k*df_entropy_stage_I_tumor$log2_pk
    df_entropy_stage_II_tumor$p_k_mult_log2_pk<-df_entropy_stage_II_tumor$p_k*df_entropy_stage_II_tumor$log2_pk
    df_entropy_stage_III_tumor$p_k_mult_log2_pk<-df_entropy_stage_III_tumor$p_k*df_entropy_stage_III_tumor$log2_pk

    # Caclulate entropy value
    Entropy_stage_I_value_Carels_normal  <-abs(sum(df_entropy_stage_I_normal$p_k_mult_log2_pk))
    Entropy_stage_II_value_Carels_normal  <-abs(sum(df_entropy_stage_II_normal$p_k_mult_log2_pk))
    Entropy_stage_III_value_Carels_normal  <-abs(sum(df_entropy_stage_III_normal$p_k_mult_log2_pk))

    
    # Caclulate entropy value
    Entropy_stage_I_value_Carels_tumor  <-abs(sum(df_entropy_stage_I_tumor$p_k_mult_log2_pk))
    Entropy_stage_II_value_Carels_tumor  <-abs(sum(df_entropy_stage_II_tumor$p_k_mult_log2_pk))
    Entropy_stage_III_value_Carels_tumor  <-abs(sum(df_entropy_stage_III_tumor$p_k_mult_log2_pk))
}





round(Entropy_stage_I_value_Carels,4)
round(Entropy_stage_II_value_Carels,4)
round(Entropy_stage_III_value_Carels,4)
