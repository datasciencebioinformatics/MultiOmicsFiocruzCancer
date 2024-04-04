A manifest data for transcriptome profiling using RNA-Seq of TCGA lung data was created in GDC portal (https://portal.gdc.cancer.gov/) at 08-03-2024. A total of M cases and N sample is listed in the manifest file. Among these, the numbers of cases from primary tumor and normal solid tissue are respectively XX and YY; while the respective numbers of samples are  ZZ and KK. From all primary diagnosis, only entries from Squamous cell carcinoma, (not otherwise specified) were kept, K cases for primary tumor and normal solid tissue, respectively; whilst XX and YY samples.

gdc_manifest.2024-03-08.txt parameters:
Primary Diagnosis=Squamous cell carcinoma, NOS

Tissue or Organ of Origin=lower lobe, lung; lung, nos; middle lobe, lung; overlapping lesion of lung; upper lobe, lung

Experimental Strategy=RNA-Seq

Data Format=tsv

Program=TCGA


## A command line to download the samples from the portal cancer database
gdc-client download -m /home/felipe/Documentos/LungPortal/gdc_manifest.2024-04-04.txt
gdc-client download -m /home/felipe/Documentos/LungPortal/gdc_manifest.2024-03-08.txt

## A R scripts to create metadata from gdc files. 
#### Inputs:
Pipeline_SquamousCellCarcinoma_CreateMetadataFromGDCFiles.R
gdc_sample_sheet.2024-03-08.tsv, clinical.txt, sample.txt, exposure.txt, merged_data_patient_info.tsv

#### Script:
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateMetadataFromGDCFiles.R")

## A bash script to create a single table for all samples.
#### Output folder :

/home/felipe/Documentos/LungPortal/samples/

#### Output files : 
unstranded.rna_seq.augmented_star_gene_counts.tsv 
unstranded.rna_seq.augmented_star_gene_counts.tsv 
stranded_first.rna_seq.augmented_star_gene_counts.tsv 
stranded_second.rna_seq.augmented_star_gene_counts.tsv
stranded_second.rna_seq.augmented_star_gene_counts.tsv
tpm_unstranded.rna_seq.augmented_star_gene_counts.tsv
fpkm_unstranded.rna_seq.augmented_star_gene_counts.tsv 
header.txt 
gene_name.txt

#### Script:

/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateTableFromFiles.sh

## A R scripts to load expression files
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadExpressionData.R")

## A R scripts to load paired sample
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateMetadataPairedSamples.R")

## A R script to analise metadata and expression data of selected samples and cases, In progress
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_AnalyseExpressionData.R")

## A R script to compute DE genes of each stage
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateDEGenesPerStage.R")

## A R script to compute Ven Diagram for DE genes
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_VeenDiagrams.R")

## A R script to compute GO terms analysis
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_FunctionalAnalyse.R")

