<!-- GETTING STARTED -->
##### gdc_manifest to download the data
##### Command line for download:
gdc-client download -m /home/felipe/Documentos/LungPortal/gdc_manifest.2024-03-08.txt
#####  File path: 
/home/felipe/Documentos/LungPortal/gdc_manifest.2024-03-08.txt

#####  gdc_sample_sheet to download the data
#####  Command line for download:
#####  File path: 
/home/felipe/LungPortal/Documentos/gdc_sample_sheet.2024-03-08.tsv

#####  metadata.cohort to download the data
#####  Command line for download:
#####  File path: 
/home/felipe/LungPortal/Documentos/metadata.cohort.2024-03-08.json

#####  Two excel files were generated from the script:
https://github.com/datasciencebioinformatics/MultiOmicsFiocruzCancer/blob/main/ProcessPortalDataFromFile.R

#####  One with the association to the primary_diagnosis
##### and the second with the association to the Diagnosis
/home/felipe/LungPortal/Documentos/output/categorical_numeric_pvalues.xlsx
/home/felipe/LungPortal/Documentos/output/association_covariables_Diagnosis.xlsx

#####  Report results to:
Media.
Vacines, treatment, Sanofi.
Laboratorio
Ciencia

We analised population demographics of lung cancer patients and found greater incidence of cases among age XX, gender XX, race XX [my expertise]. Moreover, tumor stage were considered together with treatment responsiveness, and we found treatment XX for tumors in stage YY are more responsive, while tumors at stage UU are more aggressive [Nicolas/Carlyle expertise].  White race interacts with number of ciggarrets smoked, this can be because they have this gene epressed while not in the black race [to be discussed with doctors].  Finally, genes XX, YY and ZZ in the pathways TT are further considered as putative for cancer development/treatment and so forth [to be discussed with molecular biologists]. 

#####  Noteworhy:
Be meticulous in the obtantion of data
Criteria for the number of samples is valid (Carels dataase is smaller).
Be meticulous, simple and robust.

#####  Ethic and morals:
Morality : responsibility (money, fiocruz, time, etc), reprodutibility, interdiciplinarity (medical/scientific/bioinformatic authorities)
Ethic    : meticulous with the data, result well documented, re-usability of scripts.

