# This bash script scan all downloaded .tsv files and create a table containing information for all samples.
# Output folder /home/felipe/Documentos/LungPortal/samples/:
# unstranded.rna_seq.augmented_star_gene_counts.tsv          
# stranded_first.rna_seq.augmented_star_gene_counts.tsv
# stranded_second.rna_seq.augmented_star_gene_counts.tsv
# tpm_unstranded.rna_seq.augmented_star_gene_counts.tsv
# fpkm_unstranded.rna_seq.augmented_star_gene_counts.tsv
# header.txt ENSG gene id
# gene_name.txt gene symbol

# A script to obtain all expression data into a data table
input_folder= /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/

# gene_id	            #1
# gene_name	            #2  
# gene_typec	        #3
# unstranded	        #4	
# stranded_first	    #5
# stranded_second	    #6	
# tpm_unstranded	    #7	
#fpkm_unstranded	    #6
# List of files fo gene id
rm /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/all_samples.RPKM.tsv

# For each file
ls /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/*.csv  | while read file
do
	case_id=$(cat $file |grep -v "gene_id" | head -n 1 | sed 's\.,,\\g')
	echo $case_id
	echo $case_id > "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/"$case_id".RPKM.tsv"
	cat /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/03016c31-a2eb-4a04-9e98-f10e0c64fe9e..csv | grep -v "gene_id" | grep -v ,, | sed 's/,/\t/g' >> "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/"$case_id".RPKM.tsv"
done


# Finish here, to do
# Add collumn to the table
# Check order of gene_ids
paste /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/*.tsv  > /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/all_samples.RPKM.tsv
