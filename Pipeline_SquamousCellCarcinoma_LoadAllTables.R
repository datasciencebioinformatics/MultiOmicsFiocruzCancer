# Writing mtcars data
merge_interactome_gene_symbol_file = "/home/felipe/Documentos/LungPortal/samples/merge_interactome_gene_symbol"
colData_file = "/home/felipe/Documentos/LungPortal/samples/colData.tsv"
merged_data_patient_info_file ="/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv"
unstranded_data_file = "/home/felipe/Documentos/LungPortal/samples/unstranded.tsv"

# Load table
merge_interactome_gene_symbol  <-read.table(file = merge_interactome_gene_symbol_file, sep = '\t', header = TRUE,)         #
colData                        <-read.table(file =colData_file, sep = '\t', header = TRUE,)         #
merged_data_patient_info       <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,)         #
unstranded_data                <-read.table(file = unstranded_data_file, sep = '\t', header = TRUE,)         #
unstranded_rpkm  <-read.table(file = paste(output_dir,"unstranded_rpkm.tsv",sep=""), sep = '\t', header = TRUE,)         #
