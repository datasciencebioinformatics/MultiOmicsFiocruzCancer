#####################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output/"
#####################################################################################################################
# CreateMetadataFromGDCFiles
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadRPackages.R")

# A R scripts to create tables
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateMetadataFromGDCFiles.R")

## A R scripts to load expression files
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadExpressionData.R")

# rowMeans(unstranded_rpkm)>3 
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_RPKM_Normalization_V2.R")

#####################################################################################################################
sudo rm -r "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1A_trheshold_2A_threshold_3A_threshold_4A/"
sudo rm -r "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1A_trheshold_2A_threshold_3A_threshold_4B/"
sudo rm -r "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1A_trheshold_2A_threshold_3A_threshold_4C/"
sudo rm -r "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1B_trheshold_2A_threshold_3A_threshold_4A/"
sudo rm -r "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1B_trheshold_2A_threshold_3A_threshold_4B/"
sudo rm -r "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1B_trheshold_2A_threshold_3A_threshold_4C/"
sudo rm -r "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1C_trheshold_2A_threshold_3A_threshold_4A/"
sudo rm -r "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1C_trheshold_2A_threshold_3A_threshold_4B/"
sudo rm -r "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1C_trheshold_2A_threshold_3A_threshold_4C/"

mkdir "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1A_trheshold_2A_threshold_3A_threshold_4A/"
mkdir "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1A_trheshold_2A_threshold_3A_threshold_4B/"
mkdir "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1A_trheshold_2A_threshold_3A_threshold_4C/"
mkdir "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1B_trheshold_2A_threshold_3A_threshold_4A/"
mkdir "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1B_trheshold_2A_threshold_3A_threshold_4B/"
mkdir "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1B_trheshold_2A_threshold_3A_threshold_4C/"
mkdir "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1C_trheshold_2A_threshold_3A_threshold_4A/"
mkdir "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1C_trheshold_2A_threshold_3A_threshold_4B/"
mkdir "/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1C_trheshold_2A_threshold_3A_threshold_4C/"

# Configuration 
output_folder=$output_dir
sudo rm -R $output_folder
mkdir $output_folder
#####################################################################################################################
# # Pipeline
# RPKM threshold
trheshold_1A<-3
trheshold_1B<-5
trheshold_1C<-10


# # Pipeline
# Threshold Tumor/Normal        
trheshold_2A<-1.58

# Threshold per stage 
threshold_3A<-1.58


# Threshold correlation network
threshold_4A<-0.99990
threshold_4B<-0.999995
threshold_4C<-0.99999995
#####################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1A_trheshold_2A_threshold_3A_threshold_4A/"
#output_dir="/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1A_trheshold_2A_threshold_3A_threshold_4B/"
#output_dir="/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1A_trheshold_2A_threshold_3A_threshold_4C/"

#output_dir="/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1B_trheshold_2A_threshold_3A_threshold_4A/"
#output_dir="/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1B_trheshold_2A_threshold_3A_threshold_4B/"
#output_dir="/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1B_trheshold_2A_threshold_3A_threshold_4C/"


#output_dir="/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1C_trheshold_2A_threshold_3A_threshold_4A/"
#output_dir="/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1C_trheshold_2A_threshold_3A_threshold_4B/"
#output_dir="/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1C_trheshold_2A_threshold_3A_threshold_4C/"

output_folder=output_dir
#####################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output_tumor_trheshold_1A_trheshold_2A_threshold_3A_threshold_4A/"
output_folder=output_dir
RPKM_trheshold<-3
log2foldchange_trheshold<-0.5
log2fc_threshold<-0.25
upper_weight_th<-0.85
#####################################################################################################################
# CreateMetadataFromGDCFiles
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadRPackages.R")

## A R scripts for RPKM Normalization # Done
# rowMeans(unstranded_rpkm)>3 
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadAllTables.R")

## A R scripts to apply filters
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateMetadataPairedSamplesRPKM.R")

## A R scripts to apply filters
#source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_Filters.R")
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_Filters_V2.R")

## A R scripts to DE per stage
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateDEGenesPerStageMeansFromPairedV2.R")

## A R scripts to create veen driagam
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_VeenDiagramsFromPairedUp.R")

## A R scripts to calculate Calculate Shannon Entropy
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CalculateShannonEntropy.R")

## A R scripts to calculate Calculate Shannon Entropy
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CalculateShannonEntropyFromPairedUp.R")

## A R scripts to plot cytoscape network analysis
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_NetworkAnalysisPlot.R")

## A R scripts to plot cytoscape network analysis
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_Citoscape_Networtks.R")
