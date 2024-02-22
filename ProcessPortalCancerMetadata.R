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
#####################################################################################################################                                                                                                                    #
# tsv files  are saved as txt and null character '-- replace by -                                                   #
clinical_file="/home/felipe/portal_gdc_cancer_gov/clinical.txt"                                                     #
exposure_file="/home/felipe/portal_gdc_cancer_gov/exposure.txt"                                                     #
famhisto_file="/home/felipe/portal_gdc_cancer_gov/family_history.txt"                                               #
followup_file="/home/felipe/portal_gdc_cancer_gov/follow_up.txt"                                                    #
patholog_file="/home/felipe/portal_gdc_cancer_gov/pathology_detail.txt"                                             #
#####################################################################################################################
aliquot_file="/home/felipe/portal_gdc_cancer_gov/aliquot.txt"                                                       #
analyte_file="/home/felipe/portal_gdc_cancer_gov/analyte.txt"                                                       #
portion_file="/home/felipe/portal_gdc_cancer_gov/portion.txt"                                                       #
sample_file="/home/felipe/portal_gdc_cancer_gov/sample.txt"                                                         #
slide_file="/home/felipe/portal_gdc_cancer_gov/slide.txt"                                                           #                                                                                                  
#####################################################################################################################
output_dir="/home/felipe/portal_gdc_cancer_gov/output/"                                                             #
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
pheatmap_df_tissue_or_organ_of_origin_filtered<-pheatmap(df_tissue_or_organ_of_origin_filtered,cluster_rows = FALSE,number_format = "%.0f",display_numbers=TRUE,main="Number of cases per platform")

# FindClusters_resolution
png(filename=paste(output_dir,"/Pheatmap_number_of_cases.png",sep=""), width = 16, height = 16, res=600, units = "cm")
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

# Store pathologies with more than 50 samples per covariable
fifty_sample_names_variable<-names(table(merge_all$primary_diagnosis)[which(table(merge_all$primary_diagnosis)>50)])

# Remove first line
fifty_sample_names_variable<-fifty_sample_names_variable[fifty_sample_names_variable!="-"]

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
library("DescTools")

# Recreate merge all table
merge_all <- merge(merge_clinical_exposure_fam_followup, patholog_data, by = "case_id", suffixes = c(".merge_3","patholog"), all = TRUE, no.dups=TRUE)   

# Set co-variables
covariables <- as.vector(tolower(colnames(merge_all)))

# A tabel with data for analysis
stu_data = data.frame(covariable=merge_all$age_at_index,primary_diagnosis=factor(merge_all$primary_diagnosis))	

# Remove "-"
stu_data<-stu_data[stu_data$covariable!="-",]
stu_data<-stu_data[stu_data$primary_diagnosis!="-",]

# set rownames
stu_data<-na.omit(stu_data)

# Set covariable and primary_diagnosis
covariable<-stu_data$covariable
primary_diagnosis<-stu_data$primary_diagnosis

# Regression with categorical predictors
lm_test<-summary(lm(covariable ~ primary_diagnosis))

# Take the names of 
df_names_Pr_gt_t<-data.frame(Pr_gt_t=gsub("primary_diagnosis","",  names(lm_test$coefficients[,4]), ignore.case = FALSE, perl = FALSE,  fixed = FALSE, useBytes = FALSE),n=0)

# Set rownames
rownames(df_names_Pr_gt_t)<-df_names_Pr_gt_t$Pr_gt_t
##########################################################################################################################################################################################################
library(varhandle)
# Set co-variables
# Set co-variables
covariables <- as.vector(tolower(colnames(merge_all)))
covariables<-covariables[-which(covariables=="primary_diagnosis")]
covariables<-covariables[-which(covariables=="case_id")]
covariables<-covariables[-which(covariables=="case_submitter_id.clincal")]
covariables<-covariables[-which(covariables=="project_id.clincal")]
covariables<-covariables[-which(covariables=="case_submitter_id.exposure")]
covariables<-covariables[-which(covariables=="case_submitter_id.merge_2")]
covariables<-covariables[-which(covariables=="case_submitter_id")]

# Create a matrix for categorical results                                                                                                #
df_tissue_or_organ_of_origin_categorical_pvalues <- data.frame(matrix(Inf, ncol = length(c("chisq","goodmanKruskalGamma") ), nrow = 0))  #
                                                                                                                                         #
# Create a matrix for numeric results                                                                                                    #                                                                                   #
df_tissue_or_organ_of_origin_numeric_pvalues     <- data.frame(matrix(Inf, ncol = length(unique(merge_all$primary_diagnosis)), nrow = 0))#
                                                                                                                                         #
## Set colnames                                                                                                                          #
colnames(df_tissue_or_organ_of_origin_categorical_pvalues)<-c("chisq","goodmanKruskalGamma")                                             #                                                                                                #
                                                                                                                                         #
## Set colnames                                                                                                                          #
colnames(df_tissue_or_organ_of_origin_numeric_pvalues)<-unique(merge_all$primary_diagnosis)                                              #

# For each co-variable                                                                                                                   #
for (covariable in covariables)                                                                                                          #
{
	# Recreate merge all table
	merge_all <- merge(merge_clinical_exposure_fam_followup, patholog_data, by = "case_id", suffixes = c(".merge_3","patholog"), all = TRUE, no.dups=TRUE)               #

	# Check if co-varibale is numerical or categorical
	categorical_variable<-FALSE
	categorical_variable<-sum(check.numeric(v=merge_all[,covariable], na.rm=TRUE, only.integer=FALSE, exceptions=c(""), ignore.whitespace=TRUE))==0

	# Rename column for selected variable
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
		# if categorial covariable
		if(categorical_variable)
		{		
			# Create a contingency table with the needed variables.           
			stu_data_table = table(stu_data$covariable,stu_data$primary_diagnosis) 

			# primary_diagnosis must have at least 2 entries 
			if(dim(stu_data_table)[2]>1)
			{
				# applying chisq.test() function
				pvalue   <-chisq.test(stu_data_table)$p.value		
				
				# Applying GoodmanKruskalGamma 
				goodmanKruskalGamma<-GoodmanKruskalGamma(stu_data$covariable, y = stu_data$primary_diagnosis, conf.level = NA)
	
				# Create data.frame with results
				df_results<-data.frame(chisq=pvalue,goodmanKruskalGamma=goodmanKruskalGamma)
	
				# Set rownames
				rownames(df_results)<-covariable
	
				# Add results to data frame
				df_tissue_or_organ_of_origin_categorical_pvalues<-rbind(df_tissue_or_organ_of_origin_categorical_pvalues,df_results)			

				# Caterogical variable
				print(paste(covariable," : Categorial")	)				
				write.table(df_tissue_or_organ_of_origin_categorical_pvalues, paste(output_dir,"df_tissue_or_organ_of_origin_categorical_pvalues",".csv",sep=""), append = F)
			}	
		}
		else
		{
			# Store covariable and primary_diagnosis
			covariable_n=as.numeric(stu_data[,"covariable"])
			primary_diagnosis=stu_data[,"primary_diagnosis"]

			# primary_diagnosis must have at least 2 entries 
			if(length(unique(primary_diagnosis))>1)
			{													
				# Regression with categorical predictors
				lm_test<-summary(lm(covariable_n ~ primary_diagnosis))
	
				# "Pr(>|t|)" Pr_gt_t		
				names_Pr_gt_t<-data.frame(Pr_gt_t=gsub("primary_diagnosis","",  names(lm_test$coefficients[,4]), ignore.case = FALSE, perl = FALSE,  fixed = FALSE, useBytes = FALSE))
				Pr_gt_t<-data.frame(Pr_gt_t=lm_test$coefficients[,4], ignore.case = FALSE, perl = FALSE,  fixed = FALSE, useBytes = FALSE)
				df_results_2<-data.frame(Pr_gt_t=names_Pr_gt_t,value=Pr_gt_t$Pr_gt_t)
				colnames(df_results_2)[2]<-covariable

				# Merge by names
				df_names_Pr_gt_t<-merge(df_names_Pr_gt_t,df_results_2,by="Pr_gt_t",all.x = TRUE)
				
				# Caterogical variable
				print(paste(covariable," : Numeric")	)				
				write.table(df_names_Pr_gt_t, paste(output_dir,"df_tissue_or_organ_of_origin_numeric_pvalues",".csv",sep=""), append = F)

			}
		}
	}	
}
# Set rownames
rownames(df_names_Pr_gt_t)<-df_names_Pr_gt_t$Pr_gt_t
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
# Commentary : judicious choice of cancer type is expected according to goal.                                                                                                                            #
#            : primary inspection of database can helo selecting cancer types.                                                                                                                           #
#            : number and structure of samples too can give insights into cancer types.                                                                                                                  #
# filter pathologies with at least 50 sammples
##########################################################################################################################################################################################################
# Copy initial pavalues tabe
df_tissue_or_organ_of_origin_categorical_pvalues_clone<-df_tissue_or_organ_of_origin_categorical_pvalues
df_tissue_or_organ_of_origin_numerical_pvalues_clone  <-df_names_Pr_gt_t

# set rownames
rownames(df_tissue_or_organ_of_origin_numerical_pvalues_clone) <- df_names_Pr_gt_t$Pr_gt_t

# Set a threshold for the p-values
pvalue_threshold<-000.1

#replace all values in data frame equal to 30 with 0
df_tissue_or_organ_of_origin_categorical_pvalues_clone[df_tissue_or_organ_of_origin_categorical_pvalues_clone < pvalue_threshold]     <- 0
df_tissue_or_organ_of_origin_categorical_pvalues_clone[df_tissue_or_organ_of_origin_categorical_pvalues_clone >= pvalue_threshold]    <- 1
df_tissue_or_organ_of_origin_categorical_pvalues_clone[is.na(df_tissue_or_organ_of_origin_categorical_pvalues_clone)]                 <- 0.5

#replace all values in data frame equal to 30 with 0
df_tissue_or_organ_of_origin_numerical_pvalues_clone[df_tissue_or_organ_of_origin_numerical_pvalues_clone < pvalue_threshold]     <- 0
df_tissue_or_organ_of_origin_numerical_pvalues_clone[df_tissue_or_organ_of_origin_numerical_pvalues_clone >= pvalue_threshold]    <- 1
df_tissue_or_organ_of_origin_numerical_pvalues_clone[is.na(df_tissue_or_organ_of_origin_numerical_pvalues_clone)]                 <- 0.5

# Filter up for the variables that have at least fifty
df_tissue_or_organ_of_origin_numerical_pvalues_clone<-df_tissue_or_organ_of_origin_numerical_pvalues_clone[,colnames(df_tissue_or_organ_of_origin_numerical_pvalues_clone) %in% fifty_samples_least_variable]

# Filter up for the variables that have at least fifty
df_tissue_or_organ_of_origin_categorical_pvalues_clone<-df_tissue_or_organ_of_origin_categorical_pvalues_clone[rownames(df_tissue_or_organ_of_origin_categorical_pvalues_clone) %in% fifty_samples_least_variable,]

# fifty_sample_names_variable
df_tissue_or_organ_of_origin_numerical_pvalues_clone<-df_tissue_or_organ_of_origin_numerical_pvalues_clone[fifty_sample_names_variable,]

##########################################################################################################################################################################################################
df_tissue_or_organ_of_origin_numerical_pvalues_clone_2 = as.data.frame(sapply(df_tissue_or_organ_of_origin_numerical_pvalues_clone, as.numeric))
rownames(df_tissue_or_organ_of_origin_numerical_pvalues_clone_2) <- rownames(df_tissue_or_organ_of_origin_numerical_pvalues_clone)
df_tissue_or_organ_of_origin_numerical_pvalues_clone_2<-na.omit(df_tissue_or_organ_of_origin_numerical_pvalues_clone_2)

# Pvalues
pheatmap_tissue_or_organ_of_origin_categorical_pvalues_clone<-pheatmap(df_tissue_or_organ_of_origin_categorical_pvalues_clone,cluster_rows = TRUE,display_numbers=FALSE, main = "P-values of categorical values chi.test(covariable ~ primary_diagnosis)")
pheatmap_tissue_or_organ_of_origin_numerical_pvalues_clone<-pheatmap(df_tissue_or_organ_of_origin_numerical_pvalues_clone_2,cluster_rows = TRUE,display_numbers=FALSE, main = "P-values of numeric variables lm(covariable ~ primary_diagnosis) ")
##########################################################################################################################################################################################################
# Deliverable : a) normalized number of cases per cancer type (variables that can be used accorss pathologies)
#               b) raw number of cases per cancer type (variables that can be used per pathologies)
# a) normalized number of cases per cancer type (variables that can be used accorss pathologies)
png(filename=paste(output_dir,"/pheatmap_tissue_or_organ_of_origin_categorical_pvalues_clone.png",sep=""), width = 24, height = 32, res=600, units = "cm")
	pheatmap_tissue_or_organ_of_origin_categorical_pvalues_clone
dev.off()
#b) raw number of cases per cancer type (variables that can be used per pathologies)
png(filename=paste(output_dir,"/pheatmap_tissue_or_organ_of_origin_numerical_pvalues_clone.png",sep=""), width = 32, height = 50, res=600, units = "cm")
	pheatmap_tissue_or_organ_of_origin_numerical_pvalues_clone
dev.off()
# Questions : we will focus on the co-variables (use co-variables that can be used consistently accross pathologies?
##########################################################################################################################################################################################################
# The cancer types will be anaysed independently. 
# Use varibles that contains at least 10 entries per cancer_type accross all of them.
# Use cancer_types that are complete in covariables. 
covariables<-rownames(df_tissue_or_organ_of_origin_clone)
cancer_types<-colnames(df_tissue_or_organ_of_origin_clone)

# Binary dataset
df_tissue_or_organ_of_origin_binary<-df_tissue_or_organ_of_origin_clone*0

# For each covariable
for (covariable in covariables)
{
	# For each cancer type
	for (cancer_type in cancer_types)
	{
		# Number of samples
		n_samples<-df_tissue_or_organ_of_origin_clone[covariable,cancer_type]

		# If contains at least 30 samples
		if (n_samples >= 30 )
		{
			# Assert number to one
			df_tissue_or_organ_of_origin_binary[covariable,cancer_type]<-1		
		}
	}

}
##########################################################################################################################################################################################################
# All variables
all_variables<-rownames(df_tissue_or_organ_of_origin_clone)

# Store variables to be removes
remove_variables<-c()

# For each covariable
for (covariable in rownames(df_tissue_or_organ_of_origin_clone))
{
	print(covariable)
	if(sum(df_tissue_or_organ_of_origin_clone[covariable,]==0)>0)
	{
		# Add variable to vector
		remove_variables<-c(covariable,remove_variables)
	}	
}
# Selected variables
selected_variables<-all_variables[!all_variables %in% remove_variables]

# A list to store commom samples per cancer_types
list_cancer_types<-list()

# For each cancer type
for (cancer_type in cancer_types)
{
	# Print cancer_type
	print(cancer_type)

	# Take the samples for this cancer type
	samples_subset<-merge_all[merge_all$tissue_or_organ_of_origin==cancer_type,]

	# A vector to store complom samples
	common_samples<-samples_subset$case_id

	# For each covariable
	for (covariable in selected_variables)
	{
		# Print covariables
		covariables_values<-merge_all[merge_all$tissue_or_organ_of_origin==cancer_type,covariable]
		
		# Replace empty by NA
		covariables_values[grepl("-",covariables_values)]<-NA

		# Take subset of samples
		samples_subset<-samples_subset[which(!is.na(covariables_values)),]

		# Take the intersect samples_subset and common_samples
		common_samples<-intersect(common_samples,samples_subset$case_id)				
	}
	# Store cancer types
	list_cancer_types[[cancer_type]]<-common_samples
}
#############################################################################################################
# Table to store number of samples per cancer types
number_of_samples_per_cancer_types<-data.frame(cancer_type=c(),number_of_sameples=c())

# for each cancer type 
for (cancer_type in 1:length(list_cancer_types))
{
	# Store data
	number_of_samples_per_cancer_types<-rbind(data.frame(cancer_type=names(list_cancer_types)[[cancer_type]],number_of_sameples=length(list_cancer_types[[cancer_type]])),number_of_samples_per_cancer_types)		
}
#############################################################################################################
library("stringr")
library(dplyr)

# A procedure to verify the existance of paired-samples between cancer_type and control.
# Attempt 1 - The precudere will take all samples ids from cancer_type and all the ids from control.
# store_results$cancer_type are matched against merge_all$tissue_or_organ_of_origin and the intesection of ids is calculated
# this did not work...!

# A list to store paired samples per cancer_type
paired_sample_list<-list()

# Attempt 2 - The procedure will scan each cancer_types and search the case_id in merge_all table
# for each case, the paired samples will be searched on the table 
# for each paired sample
# for each cancer_types
for (cancer_type in cancer_types)
{
	# Store case ids for each cancer type
	cancer_type_case_ids=merge_all[which(merge_all$tissue_or_organ_of_origin==cancer_type),"case_id"]

	# Stopped here
	# To finish, check paired samples
	# remarks : I am interested in the co-variable sample_type.
	# I will se the value "normal" as control and the value "tumor" as case.
	# This way, the case must contains at least one normal and one tumor do be considered paired.
	# if this condition is met, it can be there will be more than one tumor and/or more than one case.
	# for this, combination will be created for case and sample.
	
	# Take the samples for the stored cases
	cases_cancer_type<-sample_data[which(sample_data$case_id %in% cancer_type_case_ids),"case_id"]	

	# Instance paired_sample_list
	paired_sample_list[[cancer_type]]<-data.frame(normal=c(),tumor=c(),case=c())

	# For each case, a verification for paired-samples
	for (case in cases_cancer_type)
	{
		# Store tissue types
		tissue_types<-sample_data[sample_data$case_id==case,"tissue_type"]

		# Store tssue data
		tissue_data<-sample_data[sample_data$case_id==case,]		
		

		# If there is more than one 
		if(length(tissue_types)>1)
		{			
			# if vector contains at least one tumor and one normal
			if(sum(tissue_types=="Tumor")>0 && sum(tissue_types=="Normal")>0)
			{
				# Take the "Tumor sample" samples
				tumor_samples<-tissue_data[tissue_data$tissue_type=="Tumor",]

				# Take the "Normal sample" samples
				normal_samples<-tissue_data[tissue_data$tissue_type=="Normal",]

				# A filter to keep only solid samples
				tumor_solid_samples<-tumor_samples[grepl("Solid",tumor_samples$composition ),]

				# For each tumor sample
				for (tumor_solid_sample_id in tumor_solid_samples$sample_id)
				{
					# for each normal sample, compile a paired samples
					for (normal_samples_id in normal_samples$sample_id)
					{
						# Contatenate 
						paired_sample_list[[cancer_type]]<-rbind(data.frame(normal=c(normal_samples_id),tumor=c(tumor_solid_sample_id),case=case),paired_sample_list[[cancer_type]])
					}
				}	
			}
		}		 		
	}	
}
# Data frame to store variables
df_paired_samples_info=data.frame(paired_sample=c(),number_of_cases=c(),number_of_samples=c())

# for each tissue
for (paired_sample in names(paired_sample_list))
{	
	# Take variables
	df_paired_sample<-distinct(paired_sample_list[[paired_sample]])
	paired_sample_list[[paired_sample]]<-df_paired_sample
	number_of_samples=0
	number_of_cases=0
	if (dim(df_paired_sample)[1]>0)
	{
		number_of_samples=dim(df_paired_sample)[1]
		number_of_cases=length(unique(paired_sample_list[[paired_sample]][,"case"]))
	}	
	df_paired_samples_info<-rbind(df_paired_samples_info,data.frame(paired_sample=paired_sample,number_of_cases=number_of_cases,number_of_samples=number_of_samples))
}
#########################################################################################################################
# File to case pairs
case_pairs_file="/home/felipe/portal_gdc_cancer_gov/output/paired_sample_list.xlsx"
cases_file="/home/felipe/portal_gdc_cancer_gov/output/cases_sample_list.xlsx"

#Delete file if it exists
file.remove(case_pairs_file)
file.remove(cases_file)

# Write paired-cases to file
for (tissue in names(paired_sample_list))
{
	# If there is at least one elemment
	if(dim(paired_sample_list[[tissue]])[1]>0)
	{
		# Write to table
		write.xlsx(paired_sample_list[[tissue]], file=case_pairs_file, sheetName = tissue, col.names = TRUE, row.names = TRUE, append = TRUE, showNA = TRUE, password = NULL)

		# Take all cases 
		write.xlsx(data.frame(case=unique(paired_sample_list[[tissue]]$case)), file=cases_file, sheetName = tissue, col.names = TRUE, row.names = TRUE, append = TRUE, showNA = TRUE, password = NULL)						
	}	
}
#########################################################################################################################
# Verification of the meaning of the abbreviation NOS.                                                                  #
# Essential variables   : Age Sex                                                                                       #
# Interaction variables : Race                                                                                          #
# Treatment, case and control                                                                                           #
#########################################################################################################################
# A table, each tab containing one cancer type. Each at least n>30 paired samples
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
df_tissue_or_organ_of_origin_filtered<-df_tissue_or_organ_of_origin_filtered[df_tissue_or_organ_of_origin_indexes$cancer_type,]*0

# For each tissue
for (tissue in rownames(df_tissue_or_organ_of_origin_filtered))
{			
	# For each experiment
	for (experiment in colnames(df_tissue_or_organ_of_origin_filtered))
	{
		# Start case
		cases<-c()		
		
		# Take the cases
		if(dim(paired_sample_list[[tissue]])[1]>0)
		{
			cases<-paired_sample_list[[tissue]][,"case"]
		}
		
		experiments_table<-experiments_descripton[experiments_descripton$case_id %in% cases,"data_type"]
		experiments<-names(table(experiments_table))
		experiments_count<-as.vector(table(experiments_table)[experiments])

		# Assert counts in the table
		df_tissue_or_organ_of_origin_filtered[tissue,experiments]<-experiments_count			
	}	
}
df_tissue_or_organ_of_origin_filtered$number_of_cases<-0
df_tissue_or_organ_of_origin_filtered$number_of_samples<-0

# For each line
for (tissue in rownames(df_tissue_or_organ_of_origin_filtered))
{
	df_tissue_or_organ_of_origin_filtered[tissue,"number_of_cases"]<-df_paired_samples_info[df_paired_samples_info$paired_sample==tissue,"number_of_cases"]
	df_tissue_or_organ_of_origin_filtered[tissue,"number_of_samples"]<-df_paired_samples_info[df_paired_samples_info$paired_sample==tissue,"number_of_samples"]
}

df_tissue_or_organ_of_origin_filtered[,-which(colnames(df_tissue_or_organ_of_origin_filtered)=="number_of_cases")]
df_tissue_or_organ_of_origin_filtered[,-which(colnames(df_tissue_or_organ_of_origin_filtered)=="number_of_cases")]

# FindClusters_resolution
pheatmap_df_tissue_or_organ_of_origin_filtered<-pheatmap(df_tissue_or_organ_of_origin_filtered,cluster_rows = FALSE,number_format = "%.0f",display_numbers=TRUE,main="Number of paired samples per platform")

# FindClusters_resolution
png(filename=paste(output_dir,"/Pheatmap_number_of_samples.png",sep=""), width = 16, height = 16, res=600, units = "cm")
	pheatmap_df_tissue_or_organ_of_origin_filtered
dev.off()
#########################################################################################################################
# Write a script to write the cases per cancer_type. Each type per excel sheet.

