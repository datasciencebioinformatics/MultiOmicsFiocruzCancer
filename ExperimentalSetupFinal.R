#########################################################
output_dir="/home/felipe/Documentos/LungPortal/output/" #
#########################################################
threshold_rpkm=20
threshold_tumor=1.0
threshold_stage=1.58
threshold_cor=0.75

# Create output folder name
output_folder=paste(output_dir,"/trheshold1_",threshold_rpkm,"_trheshold2_",threshold_tumor,"_trheshold3_",threshold_stage,"_trheshold4_",threshold_cor,"/",sep="")
output_folder=output_dir
############################	

# Create file							
paste("#!/usr/bin/env Rscript",sep="")
paste("threshold_rpkm=",threshold_rpkm,sep="")
paste("threshold_tumor=",threshold_tumor,sep="")				
paste("threshold_stage=",threshold_stage,sep="")								
paste("threshold_cor=",threshold_cor,sep="")												
paste("output_folder=--",output_folder,"--",sep="")
paste("output_dir=--",output_folder,"--",sep="")				

source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadRPackages.R")
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadAllTables.R")
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateMetadataPairedSamplesRPKM.R")
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_Filters_V2.R")
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateDEGenesPerStageMeansFromPairedV2.R")
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_VeenDiagramsFromPairedUp.R")
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CalculateShannonEntropy.R")
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CalculateShannonEntropyFromPairedUp.R")
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_NetworkAnalysisPlot.R")
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_Citoscape_Networtks.R")
