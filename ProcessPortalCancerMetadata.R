#####################################################################################################################
# A script to compile the table descriptive os the cases from the cancer database (https://portal.gdc.cancer.gov/)  #
# Entries:                                                                                                          #
# A) tar.gz file with Cases (n=44.451).                                                                             #
#	/home/felipe/portal_gdc_cancer_gov/                                                                               #
#	- clinical.tsv                                                                                                    # 
#	- exposure.tsv                                                                                                    # 
#	- family_history.tsv                                                                                              #
#	- follow_up.tsv                                                                                                   #
#	- pathology_detail.tsv                                                                                            #
# B) json file decription of case file (n=986.114)                                                                  #
#	- files.2024-01-30.json                                                                                           #
#	- files.2024-01-30.csv                                                                                            #
# Obs. A script was created to transform json to csv at :                                                           #
#      /home/felipe/portal_gdc_cancer_gov/covnert_simply_json.sh                                                    #  
#     clinical.tsv file was modified to include data_type field from files.2024-01-30.csv                           #
#####################################################################################################################
# Set up iput files                                                                                                 #
experiment_description_file="/home/felipe/portal_gdc_cancer_gov/files.2024-01-30.csv"                               #
                                                                                                                    #
# tsv files  are saved as txt and null character '-- replace by -                                                   #
clinical_file="/home/felipe/portal_gdc_cancer_gov/clinical.txt"                                                     #
exposure_file="/home/felipe/portal_gdc_cancer_gov/exposure.txt"                                                     #
famhisto_file="/home/felipe/portal_gdc_cancer_gov/family_history.txt"                                               #
followup_file="/home/felipe/portal_gdc_cancer_gov/follow_up.txt"                                                    #
patholog_file="/home/felipe/portal_gdc_cancer_gov/pathology_detail.txt"                                             #
#####################################################################################################################
# Read all metadata files                                                                                           #
clinical_data<-read.table(file = clinical_file, sep = '\t', header = TRUE,fill=TRUE)                                #
exposure_data<-read.table(file = exposure_file, sep = '\t', header = TRUE,fill=TRUE)                                #
famhisto_data<-read.table(file = famhisto_file, sep = '\t', header = TRUE,fill=TRUE)                                #
followup_data<-read.table(file = famhisto_file, sep = '\t', header = TRUE,fill=TRUE)                                #
patholog_data<-read.table(file = patholog_file, sep = '\t', header = TRUE,fill=TRUE)                                #
####################################################################################################################################################################################################
# Merge files without ducplicating collumns                                                                                                                                                        #                                                                                                                                                                                                 #
merge_clinical_exposure                 <- merge(clinical_data, exposure_data, by = "case_id", suffixes = c(".clincal",".exposure"), all = TRUE, no.dups=TRUE)                                     #
merge_clinical_exposure_fam             <- merge(merge_clinical_exposure, famhisto_data, by = "case_id", suffixes = c(".merge_1","famhisto"), all = TRUE, no.dups=TRUE)                            #
merge_clinical_exposure_fam_followup    <- merge(merge_clinical_exposure_fam, followup_data, by = "case_id", suffixes = c(".merge_2","famhisto"), all = TRUE, no.dups=TRUE)                        #
merge_all                               <- merge(merge_clinical_exposure_fam_followup, patholog_data, by = "case_id", suffixes = c(".merge_3","patholog"), all = TRUE, no.dups=TRUE)               #
####################################################################################################################################################################################################
# Table for the description of the data and experiments                                                       #
experiments_descripton<-read.table(file = experiment_description_file, sep = '\t', header = TRUE,fill=TRUE)   #
###############################################################################################################








