# A R script to create metadata from gdc files.
# Inputs:
# gdc_sample_sheet.2024-03-08.tsv
# /home/felipe/Documentos/LungPortal/clinical.txt
# /home/felipe/Documentos/LungPortal/sample.txt
# /home/felipe/Documentos/LungPortal/exposure.txt
# Output : merged_data_patient_info.tsv
###########################################################################################################################
merged_data_patient_info_file       <- "/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv"                  #
colData_file                        <- "/home/felipe/Documentos/LungPortal/samples/colData.tsv"                           #
###########################################################################################################################
merged_data_patient_info_data      <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)#
colData_data                       <-read.table(file = colData_file, sep = '\t', header = TRUE,fill=TRUE)                 #
rownames(colData)                  <-colData$patient_id                                                                   #
#######################################################################################################################################
colData<-na.omit(colData)                                                                                                             #
merged_data_patient_info_data<-merged_data_patient_info_data[which(merged_data_patient_info_data$patient_id %in% colData$patient_id),]#
#######################################################################################################################################
# Paired samples                                                                                                         
paired_sample_df<-data.frame(normal=c(),tumor=c(),case=c())                                                              
                                                                                                                         
# For each case, find the pairs                                                                                          
for (case in unique(merged_data_patient_info_data$case))                                                                 
{                                                                                                                        
    # All samples for case id = "case"                                                                                   
    case_samples<-merged_data_patient_info_data[merged_data_patient_info_data$case==case,]                               
                                                                                                                         
    # Take the tumor samples                                                                                           
    tumor_sampĺes <-case_samples[case_samples$tissue_type=="Tumor",]
    normal_sampĺes<-case_samples[case_samples$tissue_type=="Normal",]

    # if vector contains at least one tumor and one normal
    if(length(unique(normal_sampĺes$sample_id))>0 && length(unique(tumor_sampĺes$sample_id))>0)
    {
            # For each tumor sample
            for (tumor_solid_sample_id in tumor_sampĺes$patient_id)
            {
                # for each normal sample, compile a paired samples
                for (normal_samples_id in normal_sampĺes$patient_id)
                {
                    # Contatenate                     
                    paired_sample_df<-rbind(data.frame(normal=c(normal_samples_id),tumor=c(tumor_solid_sample_id),case=case),paired_sample_df)
                }
            }                
    }
}
#######################################################################################################################################
# folchange=Expr(Stage i)/Expr(Stage ii and II)
# folchange=rowMeans(unstranded_rpkm[,paired_sample_df$tumor])-rowMeans(unstranded_rpkm[,paired_sample_df$normal])
# folchange=rowMeans(unstranded_rpkm[,paired_sample_df$tumor])/rowMeans(unstranded_rpkm[,paired_sample_df$normal])
# log2change=l/og(folchange,2)	
# Log(FC) = mean(log2(Group1)) - mean(log2(Group2))
# Paired t-test, RPKM of paired tumor/normal samples
# Plot with 15208 genes.
# Log2foldchange
log2change=rowMeans(log(unstranded_rpkm[,paired_sample_df$tumor],2))-rowMeans(log(unstranded_rpkm[,paired_sample_df$normal],2))

LOG_CONSTANT=0.001
log2change=rowMeans(log(unstranded_rpkm[,paired_sample_df$tumor]+LOG_CONSTANT,2))-rowMeans(log(unstranded_rpkm[,paired_sample_df$normal]+LOG_CONSTANT,2))

# log2change data
log2change_tumor_control=na.omit(data.frame(gene=names(log2change),log2change=log2change))

# First, the log2foldchane tumor/normal samples is used
log2change_tumor_control$Category<-"insignificant"

# First, the log2foldchane tumor/normal samples is used
log2change_tumor_control$Pvalue<-1

# For each genes in the tabe
for (gene in log2change_tumor_control$gene)
{
	# Take p-value
	log2change_tumor_control[gene,"Pvalue"]<-t.test(x=as.numeric(unstranded_rpkm[gene,paired_sample_df$tumor]), y=as.numeric(unstranded_rpkm[gene,paired_sample_df$normal]), paired = TRUE, alternative = "two.sided")$p.value	
}
# FRD 
log2change_tumor_control$FDR<-p.adjust(log2change_tumor_control$Pvalue, method="BH")

# Categorize genes if log2foldchange >= threshold_tumor
log2change_tumor_control[intersect(which(log2change_tumor_control$FDR<=0.05), which(log2change_tumor_control$log2change>=threshold_tumor)),"Category"]<-paste("Tumor genes", sep="")
#######################################################################################################################################
# Create volcano plot
p1 <- ggplot(log2change_tumor_control, aes(log2change, -log(FDR,10))) +  geom_point(size = 2/5) +  theme_bw()
	
# Adding color to differentially expressed genes (DEGs)
p1 <- ggplot(log2change_tumor_control, aes(log2change, -log(FDR),color = Category)) + geom_point(size = 2/5,aes(color = Category))  +
  xlab("log2FoldChange") + 
  ylab("-log10(padj)") +
  scale_color_manual(values = c("black", "red")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_bw() + ggtitle(paste("Paired t-test, RPKM of paired tumor/normal samples\nlog2foldchange >=",threshold_tumor, " and FRD <= 0.05", sep="")) + guides(fill="none")
#######################################################################################################################################
# FindClusters_resolution                                                                                                                                                                                                   #
png(filename=paste(output_dir,"Volcano_plot_Tumor_Normal.png",sep=""), width = 16, height = 16, res=600, units = "cm")                                                                                                    #
	p1
dev.off() 
#######################################################################################################################################
