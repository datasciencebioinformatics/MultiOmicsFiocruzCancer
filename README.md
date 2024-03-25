## A script download the samples from the portal cancer database
gdc-client download -m /home/felipe/Documentos/LungPortal/gdc_manifest.2024-03-08.txt

## A R script to create metadata from gdc files. 
Inputs:

gdc_sample_sheet.2024-03-08.tsv, clinical.txt, sample.txt, exposure.txt, merged_data_patient_info.tsv

Script:
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateMetadataFromGDCFiles.R")

#### This bash script scan all downloaded .tsv files and create a table containing information for all samples.
Output folder :

/home/felipe/Documentos/LungPortal/samples/

Output files : unstranded.rna_seq.augmented_star_gene_counts.tsv 
unstranded.rna_seq.augmented_star_gene_counts.tsv 
stranded_first.rna_seq.augmented_star_gene_counts.tsv 
stranded_second.rna_seq.augmented_star_gene_counts.tsv
stranded_second.rna_seq.augmented_star_gene_counts.tsv
tpm_unstranded.rna_seq.augmented_star_gene_counts.tsv
fpkm_unstranded.rna_seq.augmented_star_gene_counts.tsv 
header.txt 
gene_name.txt

Script:

/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateTableFromFiles.sh

##### This R script will take the TSV file (metadata), unstranded.rna_seq.augmented_star_gene_counts (rna-seq count data), the name of genes and samples already processed for the primary_diagnosis=Squamous cell carcinoma, NOS
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadExpressionData.R")

##### A script to analise metadata and expression data of selected samples and cases, In progress
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_AnalyseExpressionData.R")
