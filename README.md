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

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEgXa1QhEH9CX1_upuGzFWUswM4cpxefAFfrV6Kx7cm41r3ruwLZ_K-CwV3kfrJqw6aiIqw-ohJmY4RUZCXSeSzjeQPiXiFHm-UcuEwe995U8AxQJiv7vykRN22Xt9CCXLmUCxZGvuNJRFKRlZ7lX-JbL3QtUVxsz7ub4k287klv3imNjTyUWa_WxONE_uOn/s16000/Volcano_Plot_Normal_Tumor_Stage_stage_I.png)

## A R script to compute Ven Diagram for DE genes
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_VeenDiagrams.R")

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEimixyTBHoL9OwUq68Bu5vY7FjPjshuawWQLc2WeYZ_oOU91aMm5RChh-6zVJP5daMuC0_ZhhT6TiJBmGBxTVqDFdybCqAKQ2UBZTDeOomNZRyfb_nT6cNx1Nvl-dsDj0mARnqgsTm46WxBmySrfQVp_Y540n0mw4GdPSyZEFawGL74HGnTg_wEw6mA_gpy/s320/Veen_diagrams.png)

## A R script to compute GO terms analysis
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_FunctionalAnalyse.R")

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEiRKWY7JMqpSAl0u6oP3qb7NRjhmGSymTQnwh86JyHSPPf9tyFJSLFCNebvHMxdzjfImMDko_EEAJ85eEs9qhXfyREExte21UKpxb9lymPt2v2FUTuj4Rs2LTVjc1TcgORoOCiKv53rFIOWr_IVkbSh1WCoyUO7xbZFkP16JtinheziavDDAepAAiwuj2U/s4000/Figure_Hypothalamus_3.png)

## A R script to compute Pipeline_SquamousCellCarcinoma_PPI_Magnitude_Entropy R
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_PPI_Magnitude_Entropy.R")



## A R script to compute Create files for raw counts per stage
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateRawCountsPerStage.R")
