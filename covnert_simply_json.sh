jq '[.[] | [.data_category , .cases[].case_id, .cases[].project[] ]]' /home/felipe/Downloads/files.2024-01-30.json > /home/felipe/portal_gdc_cancer_gov/files.2024-01-30.simplified.json
cat ./files.2024-01-30.json | jq '.[] | [.data_category , .cases[].case_id, .cases[].project[] ] | join(",")' | sed 's/,/\t/g' | sed 's/"//g'  > files.2024-01-30.csv

rm -f  /home/felipe/portal_gdc_cancer_gov/files.2024-01-30.filtered.csv
cat /home/felipe/portal_gdc_cancer_gov/files.2024-01-30.csv | while read line
do
	# Take the number of cols
	ncols=$(echo $line | wc | awk -F  '\t' '{print $2}')
	
	if [[ $ncols -lt 6 ]]
	then
		echo $line >> /home/felipe/portal_gdc_cancer_gov/files.2024-01-30.filtered.csv 
	fi
done
