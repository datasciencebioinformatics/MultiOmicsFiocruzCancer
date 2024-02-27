library(pheatmap)                                                                                                   #
library("dendextend")                                                                                               #
library("ggplot2")
#####################################################################################################################
# A script to compile the table descriptive os the cases from the cancer database (https://portal.gdc.cancer.gov/)  #
# Entries:                                                                                                          #
# A) clinical.cohort.2024-02-23.tar.gz tar.gz file with Cases                                                       #
#	/home/felipe/Documentos/LungSquaGDC/                                                                        #
#	- clinical.tsv                                                                                              # 
#	- exposure.tsv                                                                                              # 
#	- family_history.tsv                                                                                        #
#	- follow_up.tsv                                                                                             #
#	- pathology_detail.tsv                                                                                      #
# B) biospecimen.cohort.2024-02-23.tar.gz tar.gz file with sample                                                   #
#	/home/felipe/Documentos/LungSquaGDC/                                                                        #
#	- aliquot.tsv                                                                                               # 
#	- analyte.tsv                                                                                               # 
#	- portion.tsv                                                                                               #
#	- sample.tsv                                                                                                #
#	- slide.tsv                                                                                                 #
#####################################################################################################################
# Load files
#####################################################################################################################                                                                                                                    #
# tsv files  are saved as txt and null character '-- replace by -                                                 #
clinical_file="/home/felipe/Documentos/LungSquaGDC/clinical.txt"                                                  #
exposure_file="/home/felipe/Documentos/LungSquaGDC/exposure.txt"                                                  #
famhisto_file="/home/felipe/Documentos/LungSquaGDC/family_history.txt"                                            #
followup_file="/home/felipe/Documentos/LungSquaGDC/follow_up.txt"                                                 #
###################################################################################################################
aliquot_file="/home/felipe/Documentos/LungSquaGDC/aliquot.txt"                                                    #
analyte_file="/home/felipe/Documentos/LungSquaGDC/analyte.txt"                                                    #
portion_file="/home/felipe/Documentos/LungSquaGDC/portion.txt"                                                    #
sample_file="/home/felipe/Documentos/LungSquaGDC/sample.txt"                                                      #
slide_file="/home/felipe/Documentos/LungSquaGDC/slide.txt"                                                        #                                                                                                  
#####################################################################################################################
# Read all metadata files                                                                                           #
clinical_data<-read.table(file = clinical_file, sep = '\t', header = TRUE,fill=TRUE)                                #
exposure_data<-read.table(file = exposure_file, sep = '\t', header = TRUE,fill=TRUE)                                #
famhisto_data<-read.table(file = famhisto_file, sep = '\t', header = TRUE,fill=TRUE)                                #
followup_data<-read.table(file = famhisto_file, sep = '\t', header = TRUE,fill=TRUE)                                #
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
merge_all                               <- merge(merge_clinical_exposure_fam, followup_data, by = "case_id", suffixes = c(".merge_2","famhisto"), all = TRUE, no.dups=TRUE)                        #
#####################################################################################################################
output_dir="/home/felipe/Documentos/LungSquaGDC/output/"                                                           #
#####################################################################################################################
library("stringr")
library(dplyr)
# A procedure to verify the existance of paired-samples between cancer_type and control.
# Attempt 1 - The precudere will take all samples ids from cancer_type and all the ids from control.
# store_results$cancer_type are matched against merge_all$tissue_or_organ_of_origin and the intesection of ids is calculated
# this did not work...!

# Take all cases
all_cases<-unique(clinical_data$case_id)

# Take tumor and normal samples
tumor_samples<-unique(sample_data[which(sample_data$sample_type=="Primary Tumor"),"sample_id"])
normal_samples<-unique(sample_data[which(sample_data$sample_type=="Blood Derived Normal"),"sample_id"])

# Take intersect
intersect(tumor_samples,normal_samples)

# Take the union
unique(union(tumor_samples,normal_samples))
####################################################################################################################
table(merge_all$primary_diagnosis)


####################################################################################################################
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
merge_all  <- merge(merge_clinical_exposure_fam, followup_data, by = "case_id", suffixes = c(".merge_2","famhisto"), all = TRUE, no.dups=TRUE)                        #

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
#covariables<-covariables[-which(covariables=="case_submitter_id.clincal")]
#covariables<-covariables[-which(covariables=="project_id.clincal")]
#covariables<-covariables[-which(covariables=="case_submitter_id.exposure")]
#covariables<-covariables[-which(covariables=="case_submitter_id.merge_2")]
#covariables<-covariables[-which(covariables=="case_submitter_id")]
#covariables<-covariables[-which(covariables=="project_id.exposure")]

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
	merge_all <- merge(merge_clinical_exposure_fam, followup_data, by = "case_id", suffixes = c(".merge_2","famhisto"), all = TRUE, no.dups=TRUE)                        #

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
			
			# At least two values
			if(sum(!is.na(covariable_n))>2)
			{
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
write.xlsx(t(df_names_Pr_gt_t), file=paste(output_dir,"categorical_numeric_pvalues",".xlsx",sep=""), sheetName = "numeric", col.names = TRUE, row.names = TRUE, append = TRUE, showNA = TRUE, password = NULL)						
write.xlsx(t(df_tissue_or_organ_of_origin_categorical_pvalues), file=paste(output_dir,"categorical_numeric_pvalues",".xlsx",sep=""), sheetName = "categorical", col.names = TRUE, row.names = TRUE, append = TRUE, showNA = TRUE, password = NULL)						
write.xlsx(t(df_names_Pr_gt_t<0.001), file=paste(output_dir,"categorical_numeric_pvalues",".xlsx",sep=""), sheetName = "numeric_b", col.names = TRUE, row.names = TRUE, append = TRUE, showNA = TRUE, password = NULL)						
write.xlsx(t(df_tissue_or_organ_of_origin_categorical_pvalues)<0.001, file=paste(output_dir,"categorical_numeric_pvalues",".xlsx",sep=""), sheetName = "categorical_b", col.names = TRUE, row.names = TRUE, append = TRUE, showNA = TRUE, password = NULL)						


##########################################################################################################################################################################################################
# A to do listt for the weekend.
library("randomForest")

merge_all <- merge(merge_clinical_exposure_fam, followup_data, by = "case_id", suffixes = c(".merge_2","famhisto"), all = TRUE, no.dups=TRUE)                        #

# Percentage of complete data
complete_data_per_variable<-data.frame(Covariable=c(),completeness=c())

# For each column, convert to numeric
for (col_bio in colnames(merge_all))
{               
        # Percentage of complete data
        complete_data_per_variable<-rbind(complete_data_per_variable,data.frame(Covariable=c(col_bio),completeness=c(sum(!is.na(merge_all[,col_bio]))/length(merge_all[,col_bio])*100)))
}
rownames(complete_data_per_variable)<-complete_data_per_variable$Covariable     

# Sort completness table
complete_data_per_variable<-complete_data_per_variable[order(complete_data_per_variable$completeness),]

# Remove variables 
complete_data_per_variable<-complete_data_per_variable[complete_data_per_variable$completeness>0,]

# Sort data
complete_data_per_variable<-complete_data_per_variable[order(complete_data_per_variable$completeness),]

# Re-factor
complete_data_per_variable$Covariable<-factor(complete_data_per_variable$Covariable,levels=complete_data_per_variable$Covariable)

# Create plot
plt <- ggplot(complete_data_per_variable) +  geom_col(aes(completeness, Covariable)) + theme_bw() + geom_vline(xintercept=90, linetype='dotted', col = 'black',linewidth=2)

# FindClusters_resolution                                                              
png(filename=paste(output_dir,"/complete_data_per_variable.png",sep=""), width = 18, height = 24, res=600, units = "cm")
        plt
dev.off()

##########################################################################################################################################################################################################
# Covariable RACE
# Create data.frame
df_covariable<-table(merge_all$primary_diagnosis,merge_all$race)

# Take name of covariables
variable_names<-rownames(df_covariable)

# Merge table
df_covariable<-cbind(df_covariable,data.frame(names=variable_names))

# Rename colnames
colnames(df_covariable)<-c("Name","Covariable","Count","Names")

# Compute percentage
df_covariable$Percentage<-0

# Compute percentage 
for (name in unique(df_covariable$Names))
{        
	df_covariable[df_covariable$Names==name,"Percentage"]<-df_covariable[df_covariable$Names==name,"Count"]/sum(df_covariable[df_covariable$Names==name,"Count"])*100
}

# FindClusters_resolution                                                              
png(filename=paste(output_dir,"/Race_primary_diagnosis_barplot.png",sep=""), width = 24, height = 16, res=600, units = "cm")
        ggplot(df_covariable,aes(x=Names,y=Percentage,fill=Covariable))+ geom_bar(stat="identity") + theme_bw()+ coord_flip() + ggtitle("Race vs. Primary diganosis")
dev.off()
##########################################################################################################################################################################################################
# Covariable ajcc_pathologic_t
# Create data.frame
df_covariable<-table(merge_all$primary_diagnosis,merge_all$ajcc_pathologic_t)

# Take name of covariables
variable_names<-rownames(df_covariable)

# Merge table
df_covariable<-cbind(df_covariable,data.frame(names=variable_names))

# Rename colnames
colnames(df_covariable)<-c("Name","Covariable","Count","Names")

# Compute percentage
df_covariable$Percentage<-0

# Compute percentage 
for (name in unique(df_covariable$Names))
{        
	df_covariable[df_covariable$Names==name,"Percentage"]<-df_covariable[df_covariable$Names==name,"Count"]/sum(df_covariable[df_covariable$Names==name,"Count"])*100
}

# FindClusters_resolution                                                              
png(filename=paste(output_dir,"/ajcc_pathologic_t_primary_diagnosis_barplot.png",sep=""), width = 24, height = 16, res=600, units = "cm")
        ggplot(df_covariable,aes(x=Names,y=Percentage,fill=Covariable))+ geom_bar(stat="identity") + theme_bw()+ coord_flip() + ggtitle("ajcc_pathologic_t vs. Primary diganosis")
dev.off()
##########################################################################################################################################################################################################
# Df table numeric
df_table_numeric<-

# Take data frames
df_1<-data.frame(merge_all[,c("days_to_death","primary_diagnosis")],Variable="days_to_death")
df_2<-data.frame(merge_all[,c("days_to_last_follow_up","primary_diagnosis")],Variable="days_to_last_follow_up")
df_3<-data.frame(merge_all[,c("cigarettes_per_day","primary_diagnosis")],Variable="cigarettes_per_day")
df_4<-data.frame(merge_all[,c("pack_years_smoked","primary_diagnosis")],Variable="pack_years_smoked")

# Take data frames
colnames(df_1)[1]<-c("Count")
colnames(df_2)[1]<-c("Count")
colnames(df_3)[1]<-c("Count")
colnames(df_4)[1]<-c("Count")

# Store all counts
df_all_counts<-rbind(df_1,df_2,df_3,df_4)

# Remove NA lines
df_table_numeric<-na.omit(df_all_counts) 

# Re-level factor
df_table_numeric$Variable<-factor(df_table_numeric$Variable,levels=c("pack_years_smoked","cigarettes_per_day","days_to_last_follow_up","days_to_death"))

# Use semi-transparent fill
p<-ggplot(df_table_numeric, aes(x=Count, fill=primary_diagnosis, color=primary_diagnosis)) +  geom_histogram(position="identity", alpha=0.5,binwidth = NULL)  + theme_bw()+ ggtitle("days_to_death vs. Primary diganosis")+ facet_wrap("Variable", ncol = 2)

# FindClusters_resolution                                                              
png(filename=paste(output_dir,"/Numeric_variables_primary_diagnosis_barplot.png",sep=""), width = 24, height = 16, res=600, units = "cm")
        ggplot(df_table_numeric, aes(x=Count, fill=primary_diagnosis, color=primary_diagnosis)) +  geom_histogram(position="identity", alpha=0.5,binwidth = NULL)  + theme_bw()+ ggtitle("Numeric variables vs. Primary diganosis")+ facet_wrap("Variable", ncol = 2)
dev.off()

##########################################################################################################################################################################################################
library(readr)
# A table containing the Transcriptome Profiling of each sample
# Read the sample tsv file and for each "Transcriptome Profiling" read the file
 # Reading the contents of TSV file using read_tsv() method
gdc_sample_sheet_file<-"/home/felipe/Documentos/LungSquaGDC/gdc_sample_sheet.2020-05-04.txt"

# Read data
gdc_sample_sheet_data<-read.table(file = gdc_sample_sheet_file, sep = '\t', header = TRUE,fill=TRUE)    
#####################################################################################################################
# Set count
count=0

# Data frame with the results
df_results<-data.frame(Gene=c(),Count=c())

# Read all metadata files                                                                                           #
transcriptome_Profiling_samples<-gdc_sample_sheet_data[which(gdc_sample_sheet_data$Data.Category=="Transcriptome Profiling"),]

# For each file_id
for (file_id in gdc_sample_sheet_data$File.ID)
{
	# Sample file id
	file_sample_id<-list.files(paste("/home/felipe/Documentos/LungSquaGDC/samples/files/",file_id,"/",sep=""),pattern =".FPKM.txt.gz")

	# If at least one character
	if (!identical(file_sample_id, character(0)))
	{

		# Increment count
		count=count+1
		
		# Set path 
		file_path<-paste("/home/felipe/Documentos/LungSquaGDC/samples/files/",file_id,"/",file_sample_id,sep="")

		# Adjust file table
		file_sample_table<-read.table(gzfile(file_path))

		# Set colnames
		file_sample_table<-data.frame(Gene=file_sample_table$V1,Count=file_sample_table$V2)

		# Set rownames
		rownames(file_sample_table)<-file_sample_table$Gene

		# Rename colllumns
		colnames(file_sample_table)[2]<-file_id

		# If count equal to one
		if (count==1)
		{
			# Adjust file_sample_table
			df_results=file_sample_table
		}else		
		{
			# Merge table
			df_results<-merge(df_results,file_sample_table,by="Gene")
		}
		gc()
	}
}
# Set rownames
rownames(df_results)<-df_results$Gene

# Remove 
df_results<-df_results[,-1]

## To write a file in Mac Roman for simple use in Mac Excel 2004/8
write.csv(df_results, file = paste(output_dir,"Transcriptome_profiling_lung.csv",sep="/"))


