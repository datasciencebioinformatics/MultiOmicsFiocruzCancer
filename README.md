Results
	- Stage-specific singnatures for squamous cell lung cancer
	- Cancer progression measured by shannon entropy
	- Stage-wise network and funtional analysis

Stage-specific signatures were detetected followed by analysis of their respective sub-interactomes. For stages I, II and III respectivelly, AA, BB and CC genes were determined as potential biomarkers, but only DD, EE and FF were without overlapping, exclusive from each stage. From the stage-specific networks, shannon entropy was calculated to test whether there is more cancer information on expression of genes in the later stages compared to the first. Entropy values of CC, GG and HH for stages I, II and III shows that pehaps there is more information about cacer in genes specific of stage III. Stage-wise funcional analysis showed this, this and this.

The log2folchange of tumor genes is defined as means(log2 of expression for tumor samples) - means(log2 of expression for normal samples). The log2folchange of stage-specific genes is defined as means(log2 of expression for Stage X samples) - means(log2 of expression for Stages Y and X samples).

To construct the gene coexpression network, low variation genes are removed based on a threshold for the percentage of genes to be maintained (75%). Then a correlation matrix is ​​constructed using Spearman rank correlation, but only the upper diagonal is kept to avoid redundant edges. Finally, a correlation threshold of XX was used to maintain edges with significant coexpression. To construct sub-interactome networks, pairwise combinations of stage-specific genes are created and then filtered to keep edges overlapping the interctome.


A manifest data for transcriptome profiling using RNA-Seq of TCGA lung data was created in GDC portal (https://portal.gdc.cancer.gov/) at 2024-04-04. A total of M cases and N sample is listed in the manifest file. Among these, the numbers of cases from primary tumor and normal solid tissue are respectively XX and YY; while the respective numbers of samples are  ZZ and KK. From all primary diagnosis, only entries from Squamous cell carcinoma, (not otherwise specified) were kept, K cases for primary tumor and normal solid tissue, respectively; whilst XX and YY samples.

gdc_manifest.2024-03-08.txt parameters:

Primary Diagnosis=Squamous cell carcinoma, NOS

Tissue or Organ of Origin=lower lobe, lung; lung, nos; middle lobe, lung; overlapping lesion of lung; upper lobe, lung

Experimental Strategy=RNA-Seq

Data Format=tsv

Program=TCGA


## A command line to download the samples from the portal cancer database
gdc-client download -m /home/felipe/Documentos/LungPortal/gdc_manifest.2024-04-04.txt

## A R scripts to create metadata from gdc files. 
Data Category=Transcriptome Profiling

File Name    = *rna_seq.augmented_star_gene_counts.tsv

Sample.Type  = Primary Tumor or Solid Tissue Normal

#### Inputs:
Pipeline_SquamousCellCarcinoma_CreateMetadataFromGDCFiles.R
gdc_sample_sheet.2024-03-08.tsv, clinical.txt, sample.txt, exposure.txt, merged_data_patient_info.tsv

#### Script:
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateMetadataFromGDCFiles.R")

506 unique sample_id's

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

#### Scripts:
/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateTableFromFiles.sh

/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateTableFromFilesMeanOfDiffRPKM.sh


# Pipeline 
## A R scripts to create tables
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateMetadataFromGDCFiles.R")

## A R scripts to load expression files
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_LoadExpressionData.R")

## A R scripts for RPKM Normalization
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_RPKM_Normalization.R")

## A R scripts to process paired samples
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateMetadataPairedSamplesRPKM.R")

## A R scripts to calculate up-regulated genes from paired samples
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CalculateUregulatedfPairedSamplesRPKM.R")

## A R scripts to DE per stage
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateDEGenesPerStageMeansFromPairedUp.R")

## A R scripts to create veen driagam
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_VeenDiagramsFromPairedUp.R")

## A R scripts to calculate Calculate Shannon Entrop
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CalculateShannonEntropyFromPairedUp.R")

## A R scripts to assess Shannon Entropy
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_AssessShannonEntropy.R")

