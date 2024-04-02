# A R script to create metadata from gdc files.
# Inputs:
# gdc_sample_sheet.2024-03-08.tsv
# /home/felipe/Documentos/LungPortal/clinical.txt
# /home/felipe/Documentos/LungPortal/sample.txt
# /home/felipe/Documentos/LungPortal/exposure.txt
# Output : merged_data_patient_info.tsv
##########################################################################################################################################################################################################
library(readr)
library("xlsx")
library(ggplot2)
##########################################################################################################################################################################################################
output_dir="/home/felipe/Documentos/LungPortal/output/"                                                             #
#####################################################################################################################
