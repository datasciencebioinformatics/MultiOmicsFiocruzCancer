#####################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output_tumor_normal_log2fc_158_per_stage_1/"                                                             #

# Threshold Tumor/Normal        
log2foldchange_trheshold<-1.0

# Threshold per stage 
log2fc_threshold<-1
#####################################################################################################################
# Final pipeline 1
# CreateMetadataFromGDCFiles
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateMetadataFromGDCFiles.R")

## A R scripts to load expression files
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadExpressionData.R")

## A R scripts for RPKM Normalization # Done 
# rowMeans(unstranded_rpkm)>3 
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_RPKM_Normalization_V2.R")

## A R scripts to apply filters
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_Filters.R")

## A R scripts to DE per stage
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateDEGenesPerStageMeansFromPairedV2.R")

## A R scripts to create veen driagam
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_VeenDiagramsFromPairedUp.R")

## A R scripts to calculate Calculate Shannon Entropy
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CalculateShannonEntropy.R")
###########################################################################################################################




#####################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output_tumor_normal_log2fc_158_per_stage_158/"                                                             #

# Threshold Tumor/Normal        
log2foldchange_trheshold<-1.0

# Threshold per stage 
log2fc_threshold<-1.58
#####################################################################################################################
# Final pipeline 2
# CreateMetadataFromGDCFiles
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateMetadataFromGDCFiles.R")

## A R scripts to load expression files
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadExpressionData.R")

## A R scripts for RPKM Normalization # Done 
# rowMeans(unstranded_rpkm)>3 
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_RPKM_Normalization_V2.R")

## A R scripts to apply filters
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_Filters.R")

## A R scripts to DE per stage
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateDEGenesPerStageMeansFromPairedV2.R")

## A R scripts to create veen driagam
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_VeenDiagramsFromPairedUp.R")

## A R scripts to calculate Calculate Shannon Entropy
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CalculateShannonEntropy.R")
###########################################################################################################################




#####################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output_tumor_normal_log2fc_158_per_stage_2/"                                                             #

# Threshold Tumor/Normal        
log2foldchange_trheshold<-1.0

# Threshold per stage 
log2fc_threshold<-2.0
#####################################################################################################################
# Final pipeline 2
# CreateMetadataFromGDCFiles
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateMetadataFromGDCFiles.R")

## A R scripts to load expression files
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadExpressionData.R")

## A R scripts for RPKM Normalization # Done 
# rowMeans(unstranded_rpkm)>3 
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_RPKM_Normalization_V2.R")

## A R scripts to apply filters
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_Filters.R")

## A R scripts to DE per stage
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateDEGenesPerStageMeansFromPairedV2.R")

## A R scripts to create veen driagam
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_VeenDiagramsFromPairedUp.R")

## A R scripts to calculate Calculate Shannon Entropy
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CalculateShannonEntropy.R")
###########################################################################################################################




#####################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output_tumor_normal_log2fc_158_per_stage_3/"                                                             #

# Threshold Tumor/Normal        
log2foldchange_trheshold<-1.0

# Threshold per stage 
log2fc_threshold<-3.0
#####################################################################################################################
# Final pipeline 2
# CreateMetadataFromGDCFiles
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateMetadataFromGDCFiles.R")

## A R scripts to load expression files
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadExpressionData.R")

## A R scripts for RPKM Normalization # Done 
# rowMeans(unstranded_rpkm)>3 
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_RPKM_Normalization_V2.R")

## A R scripts to apply filters
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_Filters.R")

## A R scripts to DE per stage
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateDEGenesPerStageMeansFromPairedV2.R")

## A R scripts to create veen driagam
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_VeenDiagramsFromPairedUp.R")

## A R scripts to calculate Calculate Shannon Entropy
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CalculateShannonEntropy.R")
###########################################################################################################################









## A R scripts to assess Shannon Entrop
#source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_AssessShannonEntropy.R")

## A R scripts to assess Shannon Entrop
#source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_AssessNetwork.R")









## A R scripts to assess Shannon Entrop
#source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_AssessShannonEntropy.R")

## A R scripts to assess Shannon Entrop
#source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_AssessNetwork.R")

