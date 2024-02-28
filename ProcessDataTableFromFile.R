##########################################################################################################################################################################################################
library(readr)
# A table containing the Transcriptome Profiling of each sample
# Read the sample tsv file and for each "Transcriptome Profiling" read the file
 # Reading the contents of TSV file using read_tsv() method
gdc_sample_sheet_file<-"/home/felipe/Documentos/LungSquaGDC/gdc_sample_sheet.2020-05-04.txt"

# Read data
gdc_sample_sheet_data<-read.table(file = gdc_sample_sheet_file, sep = '\t', header = TRUE,fill=TRUE)  

# Add collumn sample_id
gdc_sample_sheet_data$sample_submitter_id<-gdc_sample_sheet_data$Sample.ID
#####################################################################################################################
output_dir="/home/felipe/Documentos/LungSquaGDC/output/"                                                           #
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
clinical_file="/home/felipe/Documentos/LungSquaGDC/clinical.txt" 
sample_file="/home/felipe/Documentos/LungSquaGDC/sample.txt"    

# Load data
clinical_data<-read.table(file = clinical_file, sep = '\t', header = TRUE,fill=TRUE)    
sample_data<-read.table(file = sample_file, sep = '\t', header = TRUE,fill=TRUE)                                    #

# Merge data
merged_sample_clinical_data<-merge(sample_data,clinical_data,by="case_id")

# Merge tables
merged_data_patient_info<-merge(merged_sample_clinical_data,gdc_sample_sheet_data,by="sample_submitter_id")

##########################################################################################################################################
# Take the colnames
covariables<-colnames(merged_data_patient_info)

# remove id, code, date, ID, state, code
covariables<-covariables[-which(grepl("id", covariables))]
covariables<-covariables[-which(grepl("ID", covariables))]
covariables<-covariables[-which(grepl("code", covariables))]
covariables<-covariables[-which(grepl("state", covariables))]
covariables<-covariables[-which(grepl("date", covariables))]
covariables<-covariables[-which(grepl("Name", covariables))]

# Remove co-variables
covariables<-covariables[-which(covariables %in% c("File.ID","File.Name","Data.Category","Data.Type", "Project.ID" ,"Case.ID", "Sample.ID", "sample_submitter_id", "case_id","submitter_id", "project_id.y", "project_id.x", "case_submitter_id.x", "sample_id", "submitter_id", "project_id.y","sample_type_id","submitter_id.1","case_id.1","case_submitter_id.y","created_datetime.1","state.1","submitter_id.2","created_datetime.2","updated_datetime.2","Sample.Type","biospecimen_anatomic_site","biospecimen_laterality","catalog_reference","composition","current_weight","days_to_sample_procurement","diagnosis_pathologically_confirmed","distance_normal_to_tumor","treatment_id","submitter_id","created_datetime","submitter_id","diagnosis_id","pathology_report_uuid","treatment_id","submitter_id"))]

merged_data_patient_info$Diagnosis<-NA
merged_data_patient_info[factor(merged_data_patient_info$Sample.Type)=="Primary Tumor","Diagnosis"]<-1
merged_data_patient_info[factor(merged_data_patient_info$Sample.Type)=="Solid Tissue Normal","Diagnosis"]<-0

# Save results
df_results<-data.frame(Estimate=c(), Std.Error=c(), z.value=c(), pvalue=c())

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
## To write a file in Mac Roman for simple use in Mac Excel 2004/8
write.csv(df_results, file = paste(output_dir,"Association_covariables_Diagnosis.csv",sep="/"))

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
#####################################################################################################################
# Subset datainfo
gdc_sample_sheet_data_info<-gdc_sample_sheet_data[gdc_sample_sheet_data$File.ID %in% colnames(df_results),]

# Set rownames
rownames(gdc_sample_sheet_data_info)<-gdc_sample_sheet_data_info$File.ID

# Get dataInfo
datInfo<-gdc_sample_sheet_data_info[,c("Data.Category", "Data.Type", "Sample.Type")]
