##########################################################################################################################################################################################################
library(readr)
library("xlsx")
library(ggplot2)
# A table containing the Transcriptome Profiling of each sample
# Read the sample tsv file and for each "Transcriptome Profiling" read the file
 # Reading the contents of TSV file using read_tsv() method
gdc_sample_sheet_file<-"/home/felipe/Documentos/LungPortal/gdc_sample_sheet.2024-03-08.tsv"

# Read data
gdc_sample_sheet_data<-read.table(file = gdc_sample_sheet_file, sep = '\t', header = TRUE,fill=TRUE)  

# Add collumn sample_id
gdc_sample_sheet_data$sample_submitter_id<-gdc_sample_sheet_data$Sample.ID
#####################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output/"                                                           #
#####################################################################################################################
#To do:
# - Merge data info with patient info
# - Associationg between Normal/Tumor and Covariables.
# - Check criterias for selectig data e.g. "Transcriptome Profiling" only. e.g. ".FPKM.txt.gz"
# - Outputs : a) a table dataInfo + patientInfo.
# -           b) A table with expression of selected patients.
# - To plan : 1) association covariable~primary_diagnosis
# -           2) association covariable~Normal/Tumor
# - To study : Literature of association co-variables ~ cancer
#####################################################################################################################
# Set path to files
clinical_file="/home/felipe/Documentos/LungPortal/clinical.txt" 
sample_file="/home/felipe/Documentos/LungPortal/sample.txt"    
exposure_file="/home/felipe/Documentos/LungPortal/exposure.txt"                                                  #

# Load data
clinical_data<-read.table(file = clinical_file, sep = '\t', header = TRUE,fill=TRUE)    
sample_data<-read.table(file = sample_file, sep = '\t', header = TRUE,fill=TRUE)                                    #
exposure_data<-read.table(file = exposure_file, sep = '\t', header = TRUE,fill=TRUE)                                #

# Merge data
merged_sample_clinical_data<-merge(sample_data,clinical_data,by="case_id")

# Merge all
merged_sample_clinical_data<-merge(merged_sample_clinical_data,exposure_data,by="case_id")

# Merge tables
merged_data_patient_info<-merge(merged_sample_clinical_data,gdc_sample_sheet_data,by="sample_submitter_id")
#####################################################################################################################
# length(unique(merged_data_patient_info$case_id)) # Number of cases
# length(unique(merged_data_patient_info$sample_id)) # Number of samples
# sum(unique(merged_data_patient_info[,c("sample_id","Sample.Type")])[,2]=="Primary Tumor") # Number of Primary Tumor
# sum(unique(merged_data_patient_info[,c("sample_id","Sample.Type")])[,2]=="Solid Tissue Normal") # Number of Solid Tissue Normal

# Filter tumor and normal samples
primary_tumor<-merged_data_patient_info[merged_data_patient_info[,c("sample_id","Sample.Type")][,2]=="Primary Tumor",]
solid_tissue<-merged_data_patient_info[merged_data_patient_info[,c("sample_id","Sample.Type")][,2]=="Solid Tissue Normal",]
merged_data_patient_info<-rbind(primary_tumor,solid_tissue)

# length(unique(primary_tumor$sample_id))
# length(unique(solid_tissue$sample_id))

# Population demographic
# table(unique(merged_data_patient_info[,c("sample_id","primary_diagnosis")])$primary_diagnosis)
merged_data_patient_info<-merged_data_patient_info[merged_data_patient_info$primary_diagnosis=="Squamous cell carcinoma, NOS",]

# table(unique(merged_data_patient_info[,c("sample_id","ethnicity")])$ethnicity)
# table(unique(merged_data_patient_info[,c("sample_id","gender")])$gender)
# table(unique(merged_data_patient_info[,c("sample_id","vital_status")])$vital_status)
# min(merged_data_patient_info[!is.na(merged_data_patient_info$age_at_index),"age_at_index"])
# max(merged_data_patient_info[!is.na(merged_data_patient_info$age_at_index),"age_at_index"])
# mean(merged_data_patient_info[!is.na(merged_data_patient_info$age_at_index),"age_at_index"])
##########################################################################################################################################
# Take all case ids
all_samples_ids<-unique(merged_data_patient_info$sample_id)

# A data frame to store if samples are paired
df_paired_samples<-data.frame(samples=c(), paired=c())

# For each case ID, check if the sample is present in both, solid_tissue and primary_tumor
# if present in both, the PAIRED=TRUE
for (samples in all_case_ids)
{
	# If case present in case
	sample_in_tumor<-(samples %in% unique(primary_tumor$sample_id))

	# If case present in control
	sample_in_normal<-(samples %in% unique(solid_tissue$sample_id))
	
	# if paired
	paired_sample<-(sample_in_normal && sample_in_tumor)

	# Add to the database
	df_paired_samples<-rbind(data.frame(samples=samples, paired=paired_sample),df_paired_samples)		
}
# Number of paired cases
paired_samples<-df_paired_samples[df_paired_samples$paired,"samples"]

# Paired cases
merged_data_patient_info<-merged_data_patient_info[merged_data_patient_info$case_id %in% paired_cases,]

# length(unique(merged_data_patient_info$case_id)) # Number of cases
# length(unique(merged_data_patient_info$sample_id))
length(unique(merged_data_patient_info$sample_id))


##########################################################################################################################################
# Take the colnames
covariables<-colnames(merged_data_patient_info)

# remove id, code, date, ID, state, code
covariables<-covariables[-which(grepl("id", covariables))]
covariables<-covariables[-which(grepl("ID", covariables))]
covariables<-covariables[-which(grepl("code", covariables))]
covariables<-covariables[-which(grepl("state", covariables))]


# Remove co-variables
covariables<-covariables[-which(covariables %in% c("File.ID","File.Name","Data.Category","Data.Type", "Project.ID" ,"Case.ID", "Sample.ID", "sample_submitter_id", "case_id","submitter_id", "project_id.y", "project_id.x", "case_submitter_id.x", "sample_id", "submitter_id", "project_id.y","sample_type_id","submitter_id.1","case_id.1","case_submitter_id.y","created_datetime.1","state.1","submitter_id.2","created_datetime.2","updated_datetime.2","Sample.Type","biospecimen_anatomic_site","biospecimen_laterality","catalog_reference","composition","current_weight","days_to_sample_procurement","diagnosis_pathologically_confirmed","distance_normal_to_tumor","treatment_id","submitter_id","created_datetime","submitter_id","diagnosis_id","pathology_report_uuid","treatment_id","submitter_id","primary_diagnosis"))]

merged_data_patient_info$Diagnosis<-NA
merged_data_patient_info[factor(merged_data_patient_info$Sample.Type)=="Primary Tumor","Diagnosis"]<-1
merged_data_patient_info[factor(merged_data_patient_info$Sample.Type)=="Solid Tissue Normal","Diagnosis"]<-0

# Save results
df_results<-data.frame(Estimate=c(), Std.Error=c(), z.value=c(), pvalue=c())
##########################################################################################################################################

# For each covariables
for (covariable in covariables)
{
	print(covariable)

	# Try glm function
	error <- try(glm.full<-glm(formula=formula(paste(" Diagnosis ~ ",covariable,sep="")),family=binomial(link='logit'), data=merged_data_patient_info, na.action=na.omit),TRUE)

	# try-error
	if (!class(error)[1] == "try-error")
	{					
		# Produce glm model			
		# Take summary
		df_summary<-summary(glm.full)$coefficients
	
		# Rename collumns
		colnames(df_summary)<-c("Estimate","Std.Error","z.value","pvalue")
	
		# Concatenate results
		df_results<-rbind(df_results,df_summary)
	}
}
# Put star for caterogories of siginificance
for (line in rownames(df_results))
{
	for (colunmn in colnames(df_results))
	{
		value<-df_results[line,colunmn]		
		# Replace values
		replaced_value<-value<-format(as.numeric(value), scientific = TRUE,digits =4)
		df_results[line,colunmn]<-replaced_value

		print(replaced_value)
		
		if(as.numeric(value)<0.06)
		{	
			# Replace values
			df_results[line,colunmn]<-paste(replaced_value,"¬",sep="")			
		}		
		if(as.numeric(value)<0.05)
		{
			# Replace values
			df_results[line,colunmn]<-paste(replaced_value,"*",sep="")			
		}				
		if(as.numeric(value)<0.01)
		{
			# Replace values
			df_results[line,colunmn]<-paste(replaced_value,"**",sep="")			
		}				
		if(as.numeric(value)<0.001)
		{
			# Replace values
			df_results[line,colunmn]<-paste(replaced_value,"***",sep="")			
		}
	}
}

## To write a file in Mac Roman for simple use in Mac Excel 2004/8
write.xlsx(df_results, file=paste(output_dir,"association_covariables_Diagnosis",".xlsx",sep=""), sheetName = "numeric", col.names = TRUE, row.names = TRUE, append = TRUE, showNA = TRUE, password = NULL)						

#####################################################################################################################
# To Do : Download files
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
	file_sample_id<-list.files(paste("/home/felipe/Documentos/LungPortal/samples/files/",file_id,"/",sep=""),pattern =".FPKM.txt.gz")

	# If at least one character
	if (!identical(file_sample_id, character(0)))
	{

		# Increment count
		count=count+1
		
		# Set path 
		file_path<-paste("/home/felipe/Documentos/LungPortal/samples/files/",file_id,"/",file_sample_id,sep="")

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
#####################################################################################################################
library("DescTools")

# A tabel with data for analysis
stu_data = data.frame(covariable=merged_data_patient_info$age_at_index,primary_diagnosis=factor(merged_data_patient_info$primary_diagnosis))	

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

# Create a matrix for categorical results                                                                                                #
df_tissue_or_organ_of_origin_categorical_pvalues <- data.frame(matrix(Inf, ncol = length(c("chisq","goodmanKruskalGamma") ), nrow = 0))  #
                                                                                                                                         #
# Create a matrix for numeric results                                                                                                    #                                                                                   #
df_tissue_or_organ_of_origin_numeric_pvalues     <- data.frame(matrix(Inf, ncol = length(unique(merged_data_patient_info$primary_diagnosis)), nrow = 0))#

# Merge tables
merged_data_patient_info<-merge(merged_sample_clinical_data,gdc_sample_sheet_data,by="sample_submitter_id")
                                                                                                                                         #
## Set colnames                                                                                                                          #
colnames(df_tissue_or_organ_of_origin_categorical_pvalues)<-c("chisq","goodmanKruskalGamma")                                             #                                                                                                #
                                                                                                                                         #
## Set colnames                                                                                                                          #
colnames(df_tissue_or_organ_of_origin_numeric_pvalues)<-unique(merged_data_patient_info$primary_diagnosis)                               #

# Take the colnames
covariables<-colnames(merged_data_patient_info)

# remove id, code, date, ID, state, code
covariables<-covariables[-which(grepl("id", covariables))]
covariables<-covariables[-which(grepl("ID", covariables))]
covariables<-covariables[-which(grepl("code", covariables))]
covariables<-covariables[-which(grepl("state", covariables))]

# Remove co-variables
covariables<-covariables[-which(covariables %in% c("File.ID","File.Name","Data.Category","Data.Type", "Project.ID" ,"Case.ID", "Sample.ID", "sample_submitter_id", "case_id","submitter_id", "project_id.y", "project_id.x", "case_submitter_id.x", "sample_id", "submitter_id", "project_id.y","sample_type_id","submitter_id.1","case_id.1","case_submitter_id.y","created_datetime.1","state.1","submitter_id.2","created_datetime.2","updated_datetime.2","Sample.Type","biospecimen_anatomic_site","biospecimen_laterality","catalog_reference","composition","current_weight","days_to_sample_procurement","diagnosis_pathologically_confirmed","distance_normal_to_tumor","treatment_id","submitter_id","created_datetime","submitter_id","diagnosis_id","pathology_report_uuid","treatment_id","submitter_id","primary_diagnosis"))]

# For each co-variable                                                                                                                   #
for (covariable in covariables)                                                                                                          #
{
	# Recreate merge all table
	# Merge tables
	merged_data_patient_info<-merge(merged_sample_clinical_data,gdc_sample_sheet_data,by="sample_submitter_id")

	# Check if co-varibale is numerical or categorical
	categorical_variable<-FALSE
	categorical_variable<-sum(check.numeric(v=merged_data_patient_info[,covariable], na.rm=TRUE, only.integer=FALSE, exceptions=c(""), ignore.whitespace=TRUE))==0

	# Rename column for selected variable
	colnames(merged_data_patient_info)[which(colnames(merged_data_patient_info)==covariable)]<-"covariable"

	# A tabel with data for analysis
	stu_data = data.frame(covariable=merged_data_patient_info$covariable,primary_diagnosis=merged_data_patient_info$primary_diagnosis)	

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
					df_names_Pr_gt_t<-merge(df_names_Pr_gt_t,df_results_2,by="Pr_gt_t",all.x = FALSE, ,all.y = FALSE)
					
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
# Differential categorial from numerical predictors.
# Set to inf
df_names_Pr_gt_t[is.na(df_names_Pr_gt_t)] <- Inf

# Remove tw
df_names_Pr_gt_t<-df_names_Pr_gt_t[,-c(1,2)]

# Put star for caterogories of siginificance
for (line in rownames(df_names_Pr_gt_t))
{
	for (colunmn in colnames(df_names_Pr_gt_t))
	{
		value<-df_names_Pr_gt_t[line,colunmn]		
		# Replace values
		replaced_value<-value<-format(as.numeric(value), scientific = TRUE,digits =4)
		df_names_Pr_gt_t[line,colunmn]<-replaced_value

		print(replaced_value)
		
		if(as.numeric(value)<0.06)
		{	
			# Replace values
			df_names_Pr_gt_t[line,colunmn]<-paste(replaced_value,"¬",sep="")			
		}		
		if(as.numeric(value)<0.05)
		{
			# Replace values
			df_names_Pr_gt_t[line,colunmn]<-paste(replaced_value,"*",sep="")			
		}				
		if(as.numeric(value)<0.01)
		{
			# Replace values
			df_names_Pr_gt_t[line,colunmn]<-paste(replaced_value,"**",sep="")			
		}				
		if(as.numeric(value)<0.001)
		{
			# Replace values
			df_names_Pr_gt_t[line,colunmn]<-paste(replaced_value,"***",sep="")			
		}
	}
}

# To do :
# Differential categorial from numerical predictors.
# Set to inf
df_tissue_or_organ_of_origin_categorical_pvalues[is.na(df_tissue_or_organ_of_origin_categorical_pvalues)] <- Inf


# Put star for caterogories of siginificance
for (line in rownames(df_tissue_or_organ_of_origin_categorical_pvalues))
{
	for (colunmn in colnames(df_tissue_or_organ_of_origin_categorical_pvalues))
	{
		value<-df_tissue_or_organ_of_origin_categorical_pvalues[line,colunmn]
		df_tissue_or_organ_of_origin_categorical_pvalues[line,colunmn]<-round(as.numeric(value),4)

		# Replace values
		replaced_value<-value<-format(as.numeric(value), scientific = TRUE,digits =4)

		print(replaced_value)

		
		if(as.numeric(value)<0.06)
		{
			# Replace values
			df_tissue_or_organ_of_origin_categorical_pvalues[line,colunmn]<-paste(replaced_value,"¬",sep="")			
		}		
		if(as.numeric(value)<0.05)
		{
			# Replace values
			df_tissue_or_organ_of_origin_categorical_pvalues[line,colunmn]<-paste(replaced_value,"*",sep="")			
		}				
		if(as.numeric(value)<0.01)
		{
			# Replace values
			df_tissue_or_organ_of_origin_categorical_pvalues[line,colunmn]<-paste(replaced_value,"**",sep="")			
		}				
		if(as.numeric(value)<0.001)
		{
			# Replace values
			df_tissue_or_organ_of_origin_categorical_pvalues[line,colunmn]<-paste(replaced_value,"***",sep="")			
		}
	}
}

# To do :
# chi-square
# anova
t_df_names_Pr_gt_t<-t(df_names_Pr_gt_t)
t_df_names_Pr_gt_t

write.xlsx(t(df_names_Pr_gt_t), file=paste(output_dir,"categorical_numeric_pvalues",".xlsx",sep=""), sheetName = "numeric", col.names = TRUE, row.names = TRUE, append = TRUE, showNA = TRUE, password = NULL)						
write.xlsx(t(df_tissue_or_organ_of_origin_categorical_pvalues), file=paste(output_dir,"categorical_numeric_pvalues",".xlsx",sep=""), sheetName = "categorical", col.names = TRUE, row.names = TRUE, append = TRUE, showNA = TRUE, password = NULL)						


##########################################################################################################################################################################################################
# Percentage of complete data
complete_data_per_variable<-data.frame(Covariable=c(),completeness=c())

# For each column, convert to numeric
for (col_bio in colnames(merged_data_patient_info))
{               
        # Percentage of complete data
        complete_data_per_variable<-rbind(complete_data_per_variable,data.frame(Covariable=c(col_bio),completeness=c(sum(!is.na(merged_data_patient_info[,col_bio]))/length(merged_data_patient_info[,col_bio])*100)))
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
df_covariable<-table(merged_data_patient_info$primary_diagnosis,merged_data_patient_info$race)

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
df_covariable<-table(merged_data_patient_info$primary_diagnosis,merged_data_patient_info$ajcc_pathologic_t)

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
merge_all<-merged_data_patient_info

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
