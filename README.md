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

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEhFM7z2mJa-KbtVaNSSptc6FTKCtbxJBeTtVfqzrzBJMDxADSGo248GaFVC6IEC9hKx8KtRPRUFxmGw3NE9u8RXoGEAZK4-HJnWo1eyiXSS0m9m4EY2_2u9X8MhM9hXvui3VUzTzt4rmb3HcMhxPG6F-kPiucTHBH9l1OAMYtLXDcKbpFXERCiMYhoEUoh8/s16000/Panel_genes_stage1.png)

## A R script to compute Ven Diagram for DE genes
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_VeenDiagrams.R")

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEimixyTBHoL9OwUq68Bu5vY7FjPjshuawWQLc2WeYZ_oOU91aMm5RChh-6zVJP5daMuC0_ZhhT6TiJBmGBxTVqDFdybCqAKQ2UBZTDeOomNZRyfb_nT6cNx1Nvl-dsDj0mARnqgsTm46WxBmySrfQVp_Y540n0mw4GdPSyZEFawGL74HGnTg_wEw6mA_gpy/s320/Veen_diagrams.png)

## A R script to compute GO terms analysis
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_FunctionalAnalyse.R")

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEg-mysDFmd-EqHIEeV9UKjjeWBPr0XeZkHSPaGK67xX8epdNoPaJbQS70Yyl2Mv-b5Ke9YQeMNfbuM8BDnJrs1yR5Q69ZDBkrHwFdkZnsgDV6yA4yrZw0V4g5ADW8ipVFTX_oa0kP6-UsuWyCERLl64co-zX9-6_RNxAjsmesAV-dATr52NJrBJwiZzpmWB/s3779/reducedTerms_stage_I.png)

## A R script to compute Pipeline_SquamousCellCarcinoma_PPI_Magnitude_Entropy R
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_PPI_Magnitude_Entropy.R")

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEginiCcg98ILcyMjHjJmkkFkohgB7LPTxYltJrN28c3ZCo4Y8HUKIu5OgELK1YFu8TWuhRexU3Mc4JsBeMtuKTAuSeKyi6eGdPzM3AOoCgc5VcT1i8Cm0VyLEsb6-GO6MvUmlJtCMzwNWoWgDAFGPUc32dn38fGHe1zXNS53jHVDyr5ULMiy-Ab12a3J_My/w640-h448/WhatsApp%20Image%202024-04-17%20at%2011.46.09.jpeg)

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEgQCqZ7cMW5Lm8uAuQvWxdkbD2wPR88bi65huo_KBPDPbfXKPqd3vZ5jwc3PVpoFiasG94LsnIvqIe9rYYpqrGZ9zjS5c98R59SIUgsBVdIUriIwf56Ytm-SL11kWI6UUQ8EzPYW4QWMVUlabKdTY6HP5stXOgfCB1u0qmqQy4rjTF-UHFY7fljkWNKP6TO/w640-h448/PPI_vs_log2foldchange.png)


## A R script to compute Create files for raw counts per stage
source("/home/felipe/Documentos/Fiocruz/MultiOmicsFiocruzCancer/Pipeline_SquamousCellCarcinoma_CreateRawCountsPerStage.R")
