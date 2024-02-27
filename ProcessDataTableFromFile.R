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


