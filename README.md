## A command line to download the samples from the portal cancer database
gdc-client download -m /home/felipe/Documentos/LungPortal/gdc_manifest.2024-04-04.txt

## Bash scripts:
/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateTableFromFiles.sh

/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateTableFromFilesMeanOfDiffRPKM.sh

# R Pipeline
# CreateMetadataFromGDCFiles
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadRPackages.R")

# A R scripts to create tables
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateMetadataFromGDCFiles.R")

## A R scripts to load expression files
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadExpressionData.R")

# rowMeans(unstranded_rpkm)>3 
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_RPKM_Normalization_V2.R")

## Parameters 
threshold_FDR<-0.05
threshold_tumor<-0.0
threshold_rpkm<-4
threshold_stage<-1.0
output_dir<-"/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/output_dir/"

# An R script to load all packages
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadRPackages.R")

## A R scripts to load expression files
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadAllTables.R")

## A R scripts to create tables
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateMetadataPairedSamplesRPKM.R")

## A R scripts to filter the table with RPKM threshold
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_Filters_V2.R")

## A R scripts to DE per stage
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateDEGenesPerStageMeansFromPairedV3.R")

## A R scripts to create veen driagam
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_VeenDiagramsFromPairedUp.R")

## A R scripts to assess Shannon Entropy
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CalculateShannonEntropyFromPairedUp.R")
