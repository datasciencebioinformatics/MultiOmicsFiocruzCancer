#####################################################################################################################
library(ggplot2)                                                                                                    #
library(pheatmap)                                                                                                   #
library("dendextend")                                                                                               #
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
tissue_or_organ_of_origin<-unique(clinical_data$tissue_or_organ_of_origin)[unique(clinical_data$tissue_or_organ_of_origin)!="-"]                                                                                                 #
tissue_or_organ_of_origin<-tissue_or_organ_of_origin[tissue_or_organ_of_origin!=""]                                                                                                                                              #
                                                                                                                                                                                                         #
# Create a matrix                                                                                                                                                                                        #
df_tissue_or_organ_of_origin <- data.frame(matrix(0, nrow = length(tissue_or_organ_of_origin), ncol = length(data_types)))                                                                                               #
                                                                                                                                                                                                         #
# Set rownames                                                                                                                                                                                                                                                                                                                                                   #
rownames(df_tissue_or_organ_of_origin)<-tissue_or_organ_of_origin
                                                                                                                                                                                                         #
# Set colnames                                                                                                                                                                                           #
colnames(df_tissue_or_organ_of_origin)<-data_types                                                                                                                                                               #
                                                                                                                                                                                                         #
# Fill in the table                                                                                                                                                                                      #
# for each case_ids                                                                                                                                                                                      #
for (diagnosis in tissue_or_organ_of_origin)                                                                                                                                                                     #
{                                                                                                                                                                                                        #
  # Set df_cases_type_count to 1                                                                                                                                                                         #
  df_tissue_or_organ_of_origin[diagnosis,as.vector(unique(experiments_descripton[experiments_descripton$case_id %in% clinical_data$case_id[which(clinical_data$tissue_or_organ_of_origin==diagnosis)],"data_type"]))]<-1 #
}                                                                                                                                                                                                        #
##########################################################################################################################################################################################################
# Data frame to store data frame                                                                                                                                                                         #
df_diagnosis_experiment_cancerType=data.frame(diagnosis=c(),Experiment=c(),Cancer_type=c())                                                                                                              #
                                                                                                                                                                                                         #
# For each diagnosis                                                                                                                                                                                     #
for (diagnosis in colnames(df_tissue_or_organ_of_origin))                                                                                                                                                        #
{                                                                                                                                                                                                        #
  # Fill in the table                                                                                                                                                                                    #
  df_diagnosis_experiment_cancerType=rbind(data.frame(diagnosis=diagnosis,Experiment=df_tissue_or_organ_of_origin[,diagnosis],Cancer_type=rownames(df_tissue_or_organ_of_origin)),df_diagnosis_experiment_cancerType)    #
}                                                                                                                                                                                                        #
##########################################################################################################################################################################################################
my_hclust_gene <- hclust(dist(df_tissue_or_organ_of_origin), method = "complete")
df_tissue_or_organ_of_origin<-df_tissue_or_organ_of_origin[order.hclust(my_hclust_gene),]
df_tissue_or_organ_of_origin_subset1<-df_tissue_or_organ_of_origin[1:round(dim(df_tissue_or_organ_of_origin)[1]/3)-1,]
df_tissue_or_organ_of_origin_subset2<-df_tissue_or_organ_of_origin[round(dim(df_tissue_or_organ_of_origin)[1]/3):(2*round(dim(df_tissue_or_organ_of_origin)[1]/3)),]
df_tissue_or_organ_of_origin_subset3<-df_tissue_or_organ_of_origin[(2*round(dim(df_tissue_or_organ_of_origin)[1]/3)):(3*round(dim(df_tissue_or_organ_of_origin)[1]/3))-1,]

# FindClusters_resolution
pheatmap_tissue_or_organ_of_origin<-pheatmap(df_tissue_or_organ_of_origin)
pheatmap_tissue_or_organ_of_origin_subset1<-pheatmap(df_tissue_or_organ_of_origin_subset1)
pheatmap_tissue_or_organ_of_origin_subset2<-pheatmap(df_tissue_or_organ_of_origin_subset2)
pheatmap_tissue_or_organ_of_origin_subset3<-pheatmap(df_tissue_or_organ_of_origin_subset3)

# FindClusters_resolution
png(filename=paste(output_dir,"/Pheatmap_df_tissue_or_organ_of_origin.png",sep=""), width = 24, height = 36, res=600, units = "cm")
	pheatmap_tissue_or_organ_of_origin
dev.off()

# FindClusters_resolution
png(filename=paste(output_dir,"/Pheatmap_df_tissue_or_organ_of_origin_subset_1.png",sep=""), width = 24, height = 36, res=600, units = "cm")
	pheatmap_tissue_or_organ_of_origin_subset1
dev.off()

# FindClusters_resolution
png(filename=paste(output_dir,"/Pheatmap_df_tissue_or_organ_of_origin_subset_2.png",sep=""), width = 36, height = 36, res=600, units = "cm")
	pheatmap_tissue_or_organ_of_origin_subset2
dev.off()

# FindClusters_resolution
png(filename=paste(output_dir,"/Pheatmap_df_tissue_or_organ_of_origin_subset_3.png",sep=""), width = 24, height = 36, res=600, units = "cm")
	pheatmap_tissue_or_organ_of_origin_subset3
dev.off()
##########################################################################################################################################################################################################
# To DO 03-Febrtuary-2023:
#       - Filter up by:
# 	- subsets that contains Transcription Profiling, DNA Methylation and Proteome Profilling
#	- A subset of pathologies within Fiocruz interest.
#		stomach, lung, liver, kidney, breast, Thyroid, prostate
#		https://www.frontiersin.org/articles/10.3389/fgene.2019.00930/full
#               grepl("lung",rownames(df_tissue_or_organ_of_origin) )
#               rownames(df_tissue_or_organ_of_origin)[which(grepl("lung",tolower(rownames(df_tissue_or_organ_of_origin) ) ))]
#               rownames(df_tissue_or_organ_of_origin)[which(grepl("liver",tolower(rownames(df_tissue_or_organ_of_origin) ) ))]
#               rownames(df_tissue_or_organ_of_origin)[which(grepl("kidney",tolower(rownames(df_tissue_or_organ_of_origin) ) ))]
#               rownames(df_tissue_or_organ_of_origin)[which(grepl("breast",tolower(rownames(df_tissue_or_organ_of_origin) ) ))]
#               rownames(df_tissue_or_organ_of_origin)[which(grepl("thyroid",tolower(rownames(df_tissue_or_organ_of_origin) ) ))]
#               rownames(df_tissue_or_organ_of_origin)[which(grepl("prostate",tolower(rownames(df_tissue_or_organ_of_origin) ) ))]
#               rownames(df_tissue_or_organ_of_origin)[which(grepl("stomach",tolower(rownames(df_tissue_or_organ_of_origin) ) ))]
# I have downloaded the all cancer types listed data in the cancer portal. The cancer names listed in the paper Conforte et al.2019  have a different nomeclature. I have used the word describing the "tissue or organ of origin" to match names from both portal and paper.
# My initial interest in this dataset comes from identififying variables that explain linear and non-linear associations between the differentially expressed molecules, the co-variables and the diagnosis parameters.
# Take intestine cancer, for example and diet. Consider some cultures are not willing to trade the economical profit of popular diets.
##########################################################################################################################################################################################################
# A figure containing subset of pathologies from cancer portal, after filtering for pathologies in the cancer portal listed also in Conforte et al.2019.
cancer_types<-unique(c(rownames(df_tissue_or_organ_of_origin)[which(grepl("lung",tolower(rownames(df_tissue_or_organ_of_origin) ) ))],rownames(df_tissue_or_organ_of_origin)[which(grepl("liver",tolower(rownames(df_tissue_or_organ_of_origin) ) ))],rownames(df_tissue_or_organ_of_origin)[which(grepl("kidney",tolower(rownames(df_tissue_or_organ_of_origin) ) ))],rownames(df_tissue_or_organ_of_origin)[which(grepl("breast",tolower(rownames(df_tissue_or_organ_of_origin) ) ))],rownames(df_tissue_or_organ_of_origin)[which(grepl("breast",tolower(rownames(df_tissue_or_organ_of_origin) ) ))],rownames(df_tissue_or_organ_of_origin)[which(grepl("thyroid",tolower(rownames(df_tissue_or_organ_of_origin) ) ))],rownames(df_tissue_or_organ_of_origin)[which(grepl("prostate",tolower(rownames(df_tissue_or_organ_of_origin) ) ))],rownames(df_tissue_or_organ_of_origin)[which(grepl("stomach",tolower(rownames(df_tissue_or_organ_of_origin) ) ))]))

# Filter datasets
df_tissue_or_organ_of_origin_filtered<-df_tissue_or_organ_of_origin[rownames(df_tissue_or_organ_of_origin) %in% cancer_types,]
##########################################################################################################################################################################################################
# I will now compute the number of samples per tissue+experiment·
# Create table to adjust the indexes of the pheatmap
df_tissue_or_organ_of_origin_indexes<-data.frame(cancer_type=c(rownames(df_tissue_or_organ_of_origin)[which(grepl("lung",tolower(rownames(df_tissue_or_organ_of_origin) ) ))],
rownames(df_tissue_or_organ_of_origin)[which(grepl("liver",tolower(rownames(df_tissue_or_organ_of_origin) ) ))],
rownames(df_tissue_or_organ_of_origin)[which(grepl("kidney",tolower(rownames(df_tissue_or_organ_of_origin) ) ))],
rownames(df_tissue_or_organ_of_origin)[which(grepl("breast",tolower(rownames(df_tissue_or_organ_of_origin) ) ))],
rownames(df_tissue_or_organ_of_origin)[which(grepl("thyroid",tolower(rownames(df_tissue_or_organ_of_origin) ) ))],
rownames(df_tissue_or_organ_of_origin)[which(grepl("stomach",tolower(rownames(df_tissue_or_organ_of_origin) ) ))],
rownames(df_tissue_or_organ_of_origin)[which(grepl("prostate",tolower(rownames(df_tissue_or_organ_of_origin) ) ))]),index=0)

# Set indexes
df_tissue_or_organ_of_origin_indexes$index<-1:length(df_tissue_or_organ_of_origin_indexes$index)

# Re-order
df_tissue_or_organ_of_origin_filtered<-df_tissue_or_organ_of_origin_filtered[df_tissue_or_organ_of_origin_indexes$cancer_type,]

# For each tissue
for (tissue in rownames(df_tissue_or_organ_of_origin_filtered))
{			
	# For each experiment
	for (experiment in colnames(df_tissue_or_organ_of_origin_filtered))
	{
		# Take the cases and the experiments for thoses cases
		cases<-clinical_data[clinical_data$tissue_or_organ_of_origin==tissue,"case_id"]
		experiments_table<-experiments_descripton[experiments_descripton$case_id %in% cases,"data_type"]
		experiments<-names(table(experiments_table))
		experiments_count<-as.vector(table(experiments_table)[experiments])

		# Assert counts in the table
		df_tissue_or_organ_of_origin_filtered[tissue,experiments]<-experiments_count			
	}	
}

# FindClusters_resolution
pheatmap_df_tissue_or_organ_of_origin_filtered<-pheatmap(df_tissue_or_organ_of_origin_filtered,cluster_rows = FALSE,number_format = "%.0f",display_numbers=TRUE)

# FindClusters_resolution
png(filename=paste(output_dir,"/Pheatmap_number_of_samples.png",sep=""), width = 16, height = 16, res=600, units = "cm")
	pheatmap_df_tissue_or_organ_of_origin_filtered
dev.off()
##########################################################################################################################################################################################################
# A table relating all the covriables per cancer type. It will be represented the completeness of the co-variable. Depending on the visualization, data can be filtered (only co-variables with >75% completeness cn be used). Ideally, use only samples with 100% completeness.
# Take the name of all variables
all_variables<-colnames(merge_all)

# Take the name of all pathologies
all_pathologies<-rownames(df_tissue_or_organ_of_origin_filtered)

# Create a matrix                                                                                                                                                                                        #
df_tissue_or_organ_of_origin_clone <- data.frame(matrix(0, ncol = length(all_pathologies), nrow = length(all_variables)))                                                                                               #

# Set rownames                                                                                                         #
colnames(df_tissue_or_organ_of_origin_clone)<-all_pathologies                                                                                #
                                                                                                                       #
# Set colnames                                                                                                         #
rownames(df_tissue_or_organ_of_origin_clone)<-all_variables 

# For each tissue
for (tissue in all_pathologies)
{			
	# For each experiment
	for (variable in all_variables)
	{
		# Take the cases and the experiments for thoses cases
		cases<-clinical_data[clinical_data$case_id==tissue,"case_id"]

		# Take the cases and the experiments for thoses cases
		cases<-clinical_data[clinical_data$tissue_or_organ_of_origin==tissue,"case_id"]

		# Take all variables
		variables_completeness<-as.vector(merge_all[merge_all$case_id %in% cases,variable])

		# Replace empty by NA
		variables_completeness[grepl("-",variables_completeness)]<-NA

		# Count how many are different than NA
		count_variables<-sum(!is.na(variables_completeness))
		
		# Assert counts in the table
		df_tissue_or_organ_of_origin_clone[variable,tissue]<-count_variables			

	}
}	
# Filter up with the following criteria : at least 50 samples per co-variable.
df_tissue_or_organ_of_origin_clone<-df_tissue_or_organ_of_origin_clone[which(rowSums(df_tissue_or_organ_of_origin_clone)>50),]

# Store co-variables with more than 50 samples per covariable
fifty_samples_least_variable<-rownames(df_tissue_or_organ_of_origin_clone)

# Remove 
df_tissue_or_organ_of_origin_clone<-df_tissue_or_organ_of_origin_clone[-which(rownames(df_tissue_or_organ_of_origin_clone)=="case_submitter_id.clincal"),]# FindClusters_resolution

# FindClusters_resolution
pheatmap_df_tissue_or_organ_of_origin_filtered<-pheatmap(df_tissue_or_organ_of_origin_clone,cluster_rows = FALSE,number_format = "%.0f",display_numbers=TRUE, main = "raw number of cases per cancer type")

########################################################################################################################################################################################################### 
# To answer? Which samples to use from selected pathologies?
# I have interest in the slection of covaribles and also in the selection of specific cases. For the selection of cases I have no other criteria than select cases that have complete record for the variables
# of interest. For the selection of the variables of interst, I can start with all of them. However, most are incomplete. I aim at checking the completedness of the co-variables are complete for patholgies or 
# if the compledteness if overall the same for all groups. I have computed the number of cases per co-variable but the number of patiets in each cancer type is very diverse. For this I can normalize the numbers 
# per covariables taking into account the number of cases per pathology.
# I have an interest on variables that I from  my background can have insights. Read molecular biology and systems biology. And I think too is valuable for the purpose of getting insight on the specific dataset, 
# given that this dataset too is not the whole universe of data we could have for the patholgy. Read data heterogenity. But I have share with interest this with Nicolas and my post-doc ṕeer if I read something about
# cance and this variables.
# Homeworok   : read about the co-variables and cancer data type. 
# Rember      : computer science, molecular biolgy, system biology.
# Deliverable : a) normalized number of cases per cancer type (variables that can be used accorss pathologies)
#               b) raw number of cases per cancer type (variables that can be used per pathologies)
##########################################################################################################################################################################################################
# Create a matrix                                                                                                                                                                                        #
df_tissue_or_organ_of_origin_norm <- data.frame(matrix(0, ncol = length(all_pathologies), nrow = length(all_variables)))                                                                                 #
                                                                                                                                                                                                         #
# Set rownames                                                                                                                                                                                           #
colnames(df_tissue_or_organ_of_origin_norm)<-all_pathologies                                                                                                                                             #
                                                                                                                                                                                                         #
# Set colnames                                                                                                                                                                                           #
rownames(df_tissue_or_organ_of_origin_norm)<-all_variables                                                                                                                                               #
                                                                                                                                                                                                         #
# Filter up variables with more than 50 samples per pathology                                                                                                                                            #
df_tissue_or_organ_of_origin_norm<-df_tissue_or_organ_of_origin_norm[fifty_samples_least_variable,]                                                                                                      #
                                                                                                                                                                                                         #                                                                                                                                                                                                         #
# Normalize the number for each co-variables per cancer type.                                                                                                                                            #
# The total number patients of each cance type is stored.                                                                                                                                                #
# For each tissue                                                                                                                                                                                        #
for (tissue in colnames(df_tissue_or_organ_of_origin_norm))                                                                                                                                              #
{                                                                                                                                                                                                        #
	# Take the cases                                                                                                                                                                                 #
        cases<-merge_all[merge_all$tissue_or_organ_of_origin==tissue,"case_id"]                                                                                                                          #
                                                                                                                                                                                                         #
        # Subsetlect unique cases for this tissue                                                                                                                                                        #
        unique_cases<-unique(merge_all[merge_all$case_id %in% cases,"case_id"])                                                                                                                          #                                                                                                                                                                                                         #
	                                                                                                                                                                                                 #
        # Take the total number of samples                                                                                                                                                               #
        total_cases_pathology<-length(unique_cases)                                                                                                                                                      #
                                                                                                                                                                                                         #
        # For each tissue                                                                                                                                                                                #
        for (covariables in rownames(df_tissue_or_organ_of_origin_norm))                                                                                                                                                     #
        {                                                                                                                                                                                                #
		# Some of the cases are duplicate in the table because of some co-variables have multiple entries.                                                                                       #
		# For example, I understood that when the treatment changes over time the patient data will be duplicated with the new information for the treatment.                                    #
		# In a first moment, I will not split the data per variables with multiple entries per patient, instead I will use the firtst occurance of that patient.                                 #
		merge_subset<-merge_all[merge_all$case_id %in% unique_cases,]                                                                                                                            #
                                                                                                                                                                                                         #
		# For each case_id	                                                                                                                                                                 #
		merge_subset.first <- merge_subset[match(unique(merge_subset$case_id), merge_subset$case_id),]                                                                                           #
                                                                                                                                                                                                         #
		# Take all variables                                                                                                                                                                     #
                variables_completeness<-as.vector(merge_subset.first[merge_subset.first$case_id %in% unique_cases,covariables])                                                                          #
                  #
                # Replace empty by NA
                variables_completeness[grepl("-",variables_completeness)]<-NA

                # Count how many are different than NA
                count_variables<-sum(!is.na(variables_completeness))

                # Compute percentage
                if(total_cases_pathology!=0)
                {
                        # Calculate total percentage
                        percentage_count<-(count_variables/total_cases_pathology)*100

                        # Update table
                        df_tissue_or_organ_of_origin_norm[covariables,tissue]<-percentage_count
                }                       
        }
} 
# Remove 
df_tissue_or_organ_of_origin_norm<-df_tissue_or_organ_of_origin_norm[-which(rownames(df_tissue_or_organ_of_origin_norm)=="case_submitter_id.clincal"),]# FindClusters_resolution

pheatmap_df_tissue_or_organ_of_origin_norm<-pheatmap(df_tissue_or_organ_of_origin_norm,cluster_rows = FALSE,display_numbers=TRUE, main = "normalized number of cases per cancer type")
##########################################################################################################################################################################################################
# Deliverable : a) normalized number of cases per cancer type (variables that can be used accorss pathologies)
#               b) raw number of cases per cancer type (variables that can be used per pathologies)
# a) normalized number of cases per cancer type (variables that can be used accorss pathologies)
png(filename=paste(output_dir,"/Pheatmap_number_of_cases_per_cancer_type_norm.png",sep=""), width = 32, height = 32, res=600, units = "cm")
	pheatmap_df_tissue_or_organ_of_origin_norm
dev.off()
#b) raw number of cases per cancer type (variables that can be used per pathologies)
png(filename=paste(output_dir,"/Pheatmap_raw Pheatmap_number_of_cases_per_cancer_type_raw.png",sep=""), width = 32, height = 32, res=600, units = "cm")
	pheatmap_df_tissue_or_organ_of_origin_filtered
dev.off()
# Questions : we will focus on the co-variables (use co-variables that can be used consistently accross pathologies?
##########################################################################################################################################################################################################
# Take a look at the co-variables as groups
# Treatment variables
##########################################################################################################################################################################################################
# Goal : which samples can be used for a given pathology?
# Criterias :
# 1) Samples that have all selected co-variables different than null.
# Question : how to select covariables to be used in the experiment?
# Answer   : Data has 294 cancer types and 44451 cases and 258 covariables.
# I need to find a test to assess association between the covariables and the primary diagnosis. The goal is to use the model to find co-variables to be further analysed. Most of the co-variables are categorial, 
# some are numerico. The outcome is also categorial. This way, I need to find a teste that asseses categorial predictors ~ categorial outcome. Chi-square tests, regression and data cience can be used to answer 
# this. I will try to find statistical testes that are simple enought to be compared with data science vizualization.
##########################################################################################################################################################################################################
# Recreate merge all table
merge_all <- merge(merge_clinical_exposure_fam_followup, patholog_data, by = "case_id", suffixes = c(".merge_3","patholog"), all = TRUE, no.dups=TRUE)   

# Set co-variables
covariables <- as.vector(tolower(colnames(merge_all)))

# Create a matrix                                                                                                                                                                                        #
df_tissue_or_organ_of_origin_pvalues <- data.frame(matrix(Inf, ncol = length(c("chisq","anova") ), nrow = length(covariables)))                                                                                 #
df_tissue_or_organ_of_origin_xsquared <- data.frame(matrix(Inf, ncol = length(c("chisq","anova") ), nrow = length(covariables)))                                                                                 #
df_tissue_or_organ_of_origin_df <- data.frame(matrix(Inf, ncol = length(c("chisq","anova") ), nrow = length(covariables)))                                                                                 #

                                                                                                                                                                                                         ## Set rownames                                                                                                                                                                                           #
colnames(df_tissue_or_organ_of_origin_pvalues)<-c("chisq","anova")                                                                                                                                             #
colnames(df_tissue_or_organ_of_origin_xsquared)<-c("chisq","anova")                                                                                                                                             #
colnames(df_tissue_or_organ_of_origin_df)<-c("chisq","anova")                                                                                                                                             #
                                                                                                                                                                                                         #
# Set colnames                                                                                                                                                                                           #
rownames(df_tissue_or_organ_of_origin_pvalues)<-covariables
rownames(df_tissue_or_organ_of_origin_xsquared)<-covariables
rownames(df_tissue_or_organ_of_origin_df)<-covariables

# For each co-variable
for (covariable in covariables)
{
	print(covariable)	

	# Recreate merge all table
	merge_all <- merge(merge_clinical_exposure_fam_followup, patholog_data, by = "case_id", suffixes = c(".merge_3","patholog"), all = TRUE, no.dups=TRUE)               #

	# Raname column for selected variable
	colnames(merge_all)[which(colnames(merge_all)==covariable)]<-"covariable"

	# A tabel with data for analysis
	stu_data = data.frame(covariable=merge_all$covariable,primary_diagnosis=merge_all$primary_diagnosis)

	# Remove "-"
	stu_data<-stu_data[stu_data$covariable!="-",]

	# set rownames
	stu_data<-na.omit(stu_data)

	# If there is at least one non-na entry
	if(dim(stu_data)[1]>2)
	{			
		# Create a contingency table with the needed variables.           
		stu_data = table(stu_data$covariable,stu_data$primary_diagnosis) 
	
		# applying chisq.test() function
		pvalue   <-chisq.test(stu_data)$p.value
		parameter<-chisq.test(stu_data)$parameter["df"]
		xsquared <-chisq.test(stu_data)$statistic

		# Sotre results
		df_tissue_or_organ_of_origin_pvalues[covariable,"chisq"]<-pvalue
		df_tissue_or_organ_of_origin_xsquared[covariable,"chisq"]<-parameter
		df_tissue_or_organ_of_origin_df[covariable,"chisq"]<-xsquared	
	}	
}
# To do : create three tables: 
# A table for all covariables vs. all the tests, to store p-values.
# A table for all covariables vs. all the tests, to store X-squared.
# A table for all covariables vs. all the tests, to store df.



# To do :
# Differential categorial from numerical predictors.

# To do :
# chi-square
# anova
##########################################################################################################################################################################################################
