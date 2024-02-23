library(pheatmap)                                                                                                   #
library("dendextend")                                                                                               #
#####################################################################################################################
# A script to compile the table descriptive os the cases from the cancer database (https://portal.gdc.cancer.gov/)  #
# Entries:                                                                                                          #
# A) clinical.cohort.2024-02-23.tar.gz tar.gz file with Cases                                                       #
#	/home/felipe/Documentos/Thyroid_gland/                                                                            #
#	- clinical.tsv                                                                                                    # 
#	- exposure.tsv                                                                                                    # 
#	- family_history.tsv                                                                                              #
#	- follow_up.tsv                                                                                                   #
#	- pathology_detail.tsv                                                                                            #
# B) biospecimen.cohort.2024-02-23.tar.gz tar.gz file with sample                                                   #
#	/home/felipe/Documentos/Thyroid_gland/                                                                            #
#	- aliquot.tsv                                                                                                     # 
#	- analyte.tsv                                                                                                     # 
#	- portion.tsv                                                                                                     #
#	- sample.tsv                                                                                                      #
#	- slide.tsv                                                                                                       #
#####################################################################################################################
# Load files
#####################################################################################################################                                                                                                                    #
# tsv files  are saved as txt and null character '-- replace by -                                                   #
clinical_file="/home/felipe/Documentos/Thyroid_gland/clinical.tsv"                                                  #
exposure_file="/home/felipe/Documentos/Thyroid_gland/exposure.tsv"                                                  #
famhisto_file="/home/felipe/Documentos/Thyroid_gland/family_history.tsv"                                            #
followup_file="/home/felipe/Documentos/Thyroid_gland/follow_up.tsv"                                                 #
patholog_file="/home/felipe/Documentos/Thyroid_gland/pathology_detail.tsv"                                          #
#####################################################################################################################
aliquot_file="/home/felipe/Documentos/Thyroid_gland/aliquot.tsv"                                                    #
analyte_file="/home/felipe/Documentos/Thyroid_gland/analyte.tsv"                                                    #
portion_file="/home/felipe/Documentos/Thyroid_gland/portion.tsv"                                                    #
sample_file="/home/felipe/Documentos/Thyroid_gland/sample.tsv"                                                      #
slide_file="/home/felipe/Documentos/Thyroid_gland/slide.tsv"                                                        #                                                                                                  
#####################################################################################################################
# Read all metadata files                                                                                           #
clinical_data<-read.table(file = clinical_file, sep = '\t', header = TRUE,fill=TRUE)                                #
exposure_data<-read.table(file = exposure_file, sep = '\t', header = TRUE,fill=TRUE)                                #
famhisto_data<-read.table(file = famhisto_file, sep = '\t', header = TRUE,fill=TRUE)                                #
followup_data<-read.table(file = famhisto_file, sep = '\t', header = TRUE,fill=TRUE)                                #
patholog_data<-read.table(file = patholog_file, sep = '\t', header = TRUE,fill=TRUE)                                #
#####################################################################################################################
# Read all metadata files                                                                                           #
aliquot_data<-read.table(file = aliquot_file, sep = '\t', header = TRUE,fill=TRUE)                                  #
analyte_data<-read.table(file = analyte_file, sep = '\t', header = TRUE,fill=TRUE)                                  #
portion_data<-read.table(file = portion_file, sep = '\t', header = TRUE,fill=TRUE)                                  #
sample_data<-read.table(file = sample_file, sep = '\t', header = TRUE,fill=TRUE)                                    #
slide_data<-read.table(file = slide_file, sep = '\t', header = TRUE,fill=TRUE)                                      #
#####################################################################################################################
# Merge files without ducplicating collumns                                                                                                                                                        #                                                                                                                                                                                                 #
merge_clinical_exposure                 <- merge(clinical_data, exposure_data, by = "case_id", suffixes = c(".clincal",".exposure"), all = TRUE, no.dups=TRUE)                                     #
merge_clinical_exposure_fam             <- merge(merge_clinical_exposure, famhisto_data, by = "case_id", suffixes = c(".merge_1","famhisto"), all = TRUE, no.dups=TRUE)                            #
merge_clinical_exposure_fam_followup    <- merge(merge_clinical_exposure_fam, followup_data, by = "case_id", suffixes = c(".merge_2","famhisto"), all = TRUE, no.dups=TRUE)                        #
merge_all                               <- merge(merge_clinical_exposure_fam_followup, patholog_data, by = "case_id", suffixes = c(".merge_3","patholog"), all = TRUE, no.dups=TRUE)               #
output_dir="/home/felipe/portal_gdc_cancer_gov/output/"                                                             #
#####################################################################################################################


