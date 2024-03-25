## A command line to download the samples from the portal cancer database
gdc-client download -m /home/felipe/Documentos/LungPortal/gdc_manifest.2024-03-08.txt

## A R script to create metadata from gdc files. 
#### Inputs:

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

## A R scripts to create metadata and expression files
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadExpressionData.R")

## A R script to analise metadata and expression data of selected samples and cases, In progress
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_AnalyseExpressionData.R")

## A R script to DE genes of each stage
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_AnalyseEachStage.R")
