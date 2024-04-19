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
input_folder="/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/"

# List of files fo gene id
rm -f /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageI/all_samples.RPKM.tsv
rm -f  /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageI/cases.txt
rm -r /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageI/cases_files.txt
rm -r /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageI/genes_id.txt

# For each file
ls /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageI/*.csv  | while read file
do
	case_id=$(cat $file |grep -v "gene_id" | head -n 1 | sed 's\.,,\\g')
	echo $case_id >> "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageI/cases.txt" 
	echo -e $case_id > "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageI/"$case_id".RPKM.tsv"
	cat $file | grep -v "gene_id" | grep -v ,, | sed 's/,/\t/g' | sort -u >> "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageI/"$case_id".RPKM.tsv"
 	echo "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageI/"$case_id".RPKM.tsv" >> "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageI/cases_files.txt"
  	cat $file | grep -v "gene_id" | grep -v ,, | sed 's/,/\t/g' | sort -u | awk '{print $1}' | sort -u > "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageI/genes_id.txt"
done

# List of files fo gene id
rm -f /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageII/all_samples.RPKM.tsv
rm -f  /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageII/cases.txt
rm -r /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageII/cases_files.txt
# For each file
ls /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageII/*.csv  | while read file
do
	case_id=$(cat $file |grep -v "gene_id" | head -n 1 | sed 's\.,,\\g')
	echo $case_id >> "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageII/cases.txt" 
	echo -e $case_id > "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageII/"$case_id".RPKM.tsv"
	cat $file | grep -v "gene_id" | grep -v ,, | sed 's/,/\t/g' | sort -u  >> "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageII/"$case_id".RPKM.tsv"
 	echo "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageII/"$case_id".RPKM.tsv" >> "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageII/cases_files.txt"
  	cat $file | grep -v "gene_id" | grep -v ,, | sed 's/,/\t/g' | sort -u | awk '{print $1}' | sort -u > "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageII/genes_id.txt"
done

# List of files fo gene id
rm -f /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageIII/all_samples.RPKM.tsv
rm -f  /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageIII/cases.txt
rm -r /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageIII/cases_files.txt
# For each file
ls /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageIII/*.csv  | while read file
do
	case_id=$(cat $file |grep -v "gene_id" | head -n 1 | sed 's\.,,\\g')
	echo $case_id >> "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageIII/cases.txt" 
	echo -e $case_id > "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageIII/"$case_id".RPKM.tsv"
	cat $file | grep -v "gene_id" | grep -v ,, | sed 's/,/\t/g'| sort -u  >> "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageIII/"$case_id".RPKM.tsv"
 	echo "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageIII/"$case_id".RPKM.tsv" >> "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageIII/cases_files.txt"
  	cat $file | grep -v "gene_id" | grep -v ,, | sed 's/,/\t/g' | sort -u | awk '{print $1}' | sort -u > "/home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageIII/genes_id.txt"
done


# Number of genes per stage
cat /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageI/genes_id.txt | sort -u | wc
cat /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageII/genes_id.txt |  sort -u | wc
cat /home/felipe/Documentos/LungPortal/samples/RPKM_by_Carels/all_samples/stageIII/genes_id.txt | sort -u | wc


