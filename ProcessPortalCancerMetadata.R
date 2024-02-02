#####################################################################################################################
library(ggplot2)                                                                                                    #
library(pheatmap)                                                                                                   #
#####################################################################################################################
# A script to compile the table descriptive os the cases from the cancer database (https://portal.gdc.cancer.gov/)  #
# Entries:                                                                                                          #
# A) tar.gz file with Cases (n=44.451).                                                                             #
#	/home/felipe/portal_gdc_cancer_gov/                                                                         #
#	- clinical.tsv                                                                                              # 
#	- exposure.tsv                                                                                              # 
#	- family_history.tsv                                                                                        #
#	- follow_up.tsv                                                                                             #
#	- pathology_detail.tsv                                                                                      #
# B) json file decription of case file (n=986.114)                                                                  #
#	- files.2024-01-30.json                                                                                     #
#	- files.2024-01-30.csv                                                                                      #
# Obs. A script was created to transform json to csv at :                                                           #
#      /home/felipe/portal_gdc_cancer_gov/covnert_simply_json.sh                                                    #  
#     clinical.tsv file was modified to include data_type field from files.2024-01-30.csv                           #
#####################################################################################################################
# Set up iput files                                                                                                 #
# The experiment description file orginally downloaded from the Cancer Portal plattaform was filterd as follow      #
# Only entries with number of columns less or equal to five was kept                                                #
# This was done because for some cases, more than one sample was attibuted to the experiment (case_id)              #
# mostly common in the case of Biospecimen and Clinical entries.                                                    #
# The orginal file cotains 986114 entries while the filtered file contains 985223.                                  #
experiment_description_file="/home/felipe/portal_gdc_cancer_gov/files.2024-01-30.filtered.csv"                      #
                                                                                                                    #
# tsv files  are saved as txt and null character '-- replace by -                                                   #
clinical_file="/home/felipe/portal_gdc_cancer_gov/clinical.txt"                                                     #
exposure_file="/home/felipe/portal_gdc_cancer_gov/exposure.txt"                                                     #
famhisto_file="/home/felipe/portal_gdc_cancer_gov/family_history.txt"                                               #
followup_file="/home/felipe/portal_gdc_cancer_gov/follow_up.txt"                                                    #
patholog_file="/home/felipe/portal_gdc_cancer_gov/pathology_detail.txt"                                             #
#####################################################################################################################
output_dir="/home/felipe/portal_gdc_cancer_gov/output/"                                                             #
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
# Table for the description of the data and experiments                                                                #
experiments_descripton<-unique(read.table(file = experiment_description_file, sep = '\t', header = FALSE,fill=TRUE))   #
                                                                                                                       #
# Rename colnames                                                                                                      #
colnames(experiments_descripton)<-c("data_type","case_id","project")                                                   #
########################################################################################################################
# Table : Rows = cases_id vs. Cols = data_type                                                                         #
# Clean data_types by removing entries in this collumn that do not belong to data_types.                               #
data_types<-unique(experiments_descripton$data_type)[unique(experiments_descripton$data_type)!="CGCI-BLGSP"]           #
                                                                                                                       #
# Store case_ids                                                                                                       #
case_ids<-unique(experiments_descripton$case_id)                                                                       #
                                                                                                                       #
# Create a matrix                                                                                                      #
df_cases_type_count <- data.frame(matrix(0, nrow = length(case_ids), ncol = length(data_types)))                       #
                                                                                                                       #
# Set rownames                                                                                                         #
rownames(df_cases_type_count)<-case_ids                                                                                #
                                                                                                                       #
# Set colnames                                                                                                         #
colnames(df_cases_type_count)<-data_types                                                                              #
                                                                                                                       #
# Fill in the table                                                                                                    #
# for each case_ids                                                                                                    #
for (case_id in case_ids)                                                                                              #
{                                                                                                                      #
  # Store data_types_per_case                                                                                          #
  data_types_per_case<-unique(experiments_descripton[experiments_descripton$case_id==case_id,"data_type"])             #
                                                                                                                       #
  # Clean data_types_per_case by removing entries in this collumn that do not belong to data_types.                    #
  data_types_per_case<-data_types_per_case[data_types_per_case!="CGCI-BLGSP"]                                          #
                                                                                                                       #
  # Set df_cases_type_count to 1                                                                                       #
  df_cases_type_count[case_id,data_types_per_case]<-1                                                                  #
}                                                                                                                      #
########################################################################################################################################################################################################## 
primary_diagnosis<-unique(clinical_data$primary_diagnosis)[unique(clinical_data$primary_diagnosis)!="-"]                                                                                                 #
primary_diagnosis<-primary_diagnosis[primary_diagnosis!=""]                                                                                                                                              #
                                                                                                                                                                                                         #
# Create a matrix                                                                                                                                                                                        #
df_primary_diagnosis <- data.frame(matrix(0, nrow = length(primary_diagnosis), ncol = length(data_types)))                                                                                               #
                                                                                                                                                                                                         #
# Set rownames                                                                                                                                                                                           #
rownames(df_primary_diagnosis)<-primary_diagnosis                                                                                                                                                        #
                                                                                                                                                                                                         #
# Set colnames                                                                                                                                                                                           #
colnames(df_primary_diagnosis)<-data_types                                                                                                                                                               #
                                                                                                                                                                                                         #
# Fill in the table                                                                                                                                                                                      #
# for each case_ids                                                                                                                                                                                      #
for (diagnosis in primary_diagnosis)                                                                                                                                                                     #
{                                                                                                                                                                                                        #
  # Set df_cases_type_count to 1                                                                                                                                                                         #
  df_primary_diagnosis[diagnosis,as.vector(unique(experiments_descripton[experiments_descripton$case_id %in% clinical_data$case_id[which(clinical_data$primary_diagnosis==diagnosis)],"data_type"]))]<-1 #
}                                                                                                                                                                                                        #
##########################################################################################################################################################################################################
# Data frame to store data frame                                                                                                                                                                         #
df_diagnosis_experiment_cancerType=data.frame(diagnosis=c(),Experiment=c(),Cancer_type=c())                                                                                                              #
                                                                                                                                                                                                         #
# For each diagnosis                                                                                                                                                                                     #
for (diagnosis in colnames(df_primary_diagnosis))                                                                                                                                                        #
{                                                                                                                                                                                                        #
  # Fill in the table                                                                                                                                                                                    #
  df_diagnosis_experiment_cancerType=rbind(data.frame(diagnosis=diagnosis,Experiment=df_primary_diagnosis[,diagnosis],Cancer_type=rownames(df_primary_diagnosis)),df_diagnosis_experiment_cancerType)    #
}                                                                                                                                                                                                        #
##########################################################################################################################################################################################################
# FindClusters_resolution
png(filename=paste(output_dir,"Heatmap_primary_diagnosis.png",sep=""), width = 48, height = 64, res=600, units = "cm")
	pheatmap(df_primary_diagnosis)
dev.off()
##########################################################################################################################################################################################################
