# A script to obtain all expression data into a data table
input_folder= /home/felipe/Documentos/LungPortal/samples/


gene_id	            #1
gene_name	        #2  
gene_typec	        #3
unstranded	        #4	
stranded_first	    #5
stranded_second	    #6	
tpm_unstranded	    #7	
fpkm_unstranded	    #6

# List of files fo gene id

# For each file
ls /home/felipe/Documentos/LungPortal/samples/*/*.rna_seq.augmented_star_gene_counts.tsv  | while read file
do
    folder_name=$(echo $file  | sed 's/\// /g' | awk '{print "/"$1"/"$2"/"$3"/"$4"/"$5"/"$6"/"}')
    file_name=$(echo $file  | sed 's/\// /g' | awk '{print $7}')gene_id.txt
    file_path=$(echo $folder_name"/"$file_name)
    echo $file_path
    cat $file_path | grep -v GENCODE | grep -v gene_id | grep -v N_unmapped | grep -v N_multimapping | grep -v N_noFeature | grep -v N_ambiguous | awk '{print $1}' > $folder_name$file_name".gene_id.txt"
    cat $file_path | grep -v GENCODE | grep -v gene_id | grep -v N_unmapped | grep -v N_multimapping | grep -v N_noFeature | grep -v N_ambiguous | awk '{print $2}' > $folder_name$file_name".gene_name.txt"
    cat $file_path | grep -v GENCODE | grep -v gene_id | grep -v N_unmapped | grep -v N_multimapping | grep -v N_noFeature | grep -v N_ambiguous | awk '{print $3}' > $folder_name$file_name".gene_typec.txt"
    cat $file_path | grep -v GENCODE | grep -v gene_id | grep -v N_unmapped | grep -v N_multimapping | grep -v N_noFeature | grep -v N_ambiguous | awk '{print $4}' > $folder_name$file_name".unstranded.txt"
    cat $file_path | grep -v GENCODE | grep -v gene_id | grep -v N_unmapped | grep -v N_multimapping | grep -v N_noFeature | grep -v N_ambiguous | awk '{print $5}' > $folder_name$file_name".stranded_first.txt"
    cat $file_path | grep -v GENCODE | grep -v gene_id | grep -v N_unmapped | grep -v N_multimapping | grep -v N_noFeature | grep -v N_ambiguous | awk '{print $6}' > $folder_name$file_name".stranded_second.txt"
    cat $file_path | grep -v GENCODE | grep -v gene_id | grep -v N_unmapped | grep -v N_multimapping | grep -v N_noFeature | grep -v N_ambiguous | awk '{print $7}' > $folder_name$file_name".tpm_unstranded.txt"
    cat $file_path | grep -v GENCODE | grep -v gene_id | grep -v N_unmapped | grep -v N_multimapping | grep -v N_noFeature | grep -v N_ambiguous | awk '{print $8}' > $folder_name$file_name".fpkm_unstranded.txt"

done
# Finish here, to do
# Add collumn to the table
# I stopped here 

# Paste all samples
# paste /home/felipe/Documentos/LungPortal/samples/*/*."gene_id.txt"
# paste /home/felipe/Documentos/LungPortal/samples/*/*."gene_name.txt"
# paste /home/felipe/Documentos/LungPortal/samples/*/*."gene_typec.txt"
paste /home/felipe/Documentos/LungPortal/samples/*/*."unstranded.txt"      > /home/felipe/Documentos/LungPortal/samples/unstranded.rna_seq.augmented_star_gene_counts.tmp
paste /home/felipe/Documentos/LungPortal/samples/*/*."stranded_first.txt"  > /home/felipe/Documentos/LungPortal/samples/stranded_first.rna_seq.augmented_star_gene_counts.tmp
paste /home/felipe/Documentos/LungPortal/samples/*/*."tpm_unstranded.txt"  > /home/felipe/Documentos/LungPortal/samples/tpm_unstranded.rna_seq.augmented_star_gene_counts.tmp
paste /home/felipe/Documentos/LungPortal/samples/*/*."fpkm_unstranded.txt" > /home/felipe/Documentos/LungPortal/samples/fpkm_unstranded.rna_seq.augmented_star_gene_counts.tmp

# Add gene_id to file
paste /home/felipe/Documentos/LungPortal/samples/*/*."gene_id.txt" /home/felipe/Documentos/LungPortal/samples/unstranded.rna_seq.augmented_star_gene_counts.tmp      > /home/felipe/Documentos/LungPortal/samples/unstranded.rna_seq.augmented_star_gene_counts.tsv
paste /home/felipe/Documentos/LungPortal/samples/*/*."gene_id.txt" /home/felipe/Documentos/LungPortal/samples/stranded_first.rna_seq.augmented_star_gene_counts.tmp  > /home/felipe/Documentos/LungPortal/samples/stranded_first.rna_seq.augmented_star_gene_counts.tsv
paste /home/felipe/Documentos/LungPortal/samples/*/*."gene_id.txt" /home/felipe/Documentos/LungPortal/samples/tpm_unstranded.rna_seq.augmented_star_gene_counts.tmp  > /home/felipe/Documentos/LungPortal/samples/tpm_unstranded.rna_seq.augmented_star_gene_counts.tsv
paste /home/felipe/Documentos/LungPortal/samples/*/*."gene_id.txt" /home/felipe/Documentos/LungPortal/samples/fpkm_unstranded.rna_seq.augmented_star_gene_counts.tmp > /home/felipe/Documentos/LungPortal/samples/fpkm_unstranded.rna_seq.augmented_star_gene_counts.tsv
