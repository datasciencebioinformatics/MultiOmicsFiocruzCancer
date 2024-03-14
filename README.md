<!-- GETTING STARTED -->
### Portal Cancer download criteria (create gdc_manifest)
EXPERIMENTAL_STRATEGY  = RNA-Seq

TISSUE_OR_ORGAN_ORIGIN = lower lobe, lung / lung, NOS / middle lobe, lung / upper lobe, lung / overlapping lesion of lungs

DATA_CATEGORY          = transcriptome profiling

WORKFLOW_TYPE          = STAR - Counts

DATA_FORMAT            = TSV

PROGRAM                = TCGA

### Download command line (cases in gdc_manifest)
gdc-client download -m /home/felipe/Documentos/LungPortal/gdc_manifest.2024-03-08.txt



### Filter data (primary_diagnosis==Squamous cell carcinoma, NOS)
Squamous cell carcinoma, NOS  samples (238 cases, 570 samples = 478 tumor + 92 normal), LUAD database



###  Tables
There are 238 cases, 570 samples, 478 tumor, 92 normal

/home/felipe/Documentos/LungPortal/samples/unstranded.rna_seq.gene_counts.tsv

/home/felipe/Documentos/LungPortal/samples/patient.metadata.tsv





