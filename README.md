Lassocox and Cross-validation using the Lassocox package
This repository contains code for performing lasso Cox regression analysis in R using normalized data. The analysis involves two main inputs:

## 1. Gene Expression Matrix 
The genotype expression matrix structure consists of the following format:

First column: IDs (Ensembl, gene, etc.)
Subsequent columns: Each column represents an individual sample with their corresponding gene expression levels.
Example:

```
	Ensembl_ID	TCGA1		TCGA2		TCGA3		TCGA4	
ENSG00000000003.13	11.04916787	9.665335917	10.01262454	11.14784089	
ENSG00000000460.15	7.247927513	6.426264755	7.930737338	7.247927513	
ENSG00000000971.14	6.918863237	6.50779464	7.499845887	7.977279923	
ENSG00000001036.12	10.666224	9.826548487	12.07881795	11.07347215	
ENSG00000001084.9	9.30833903	9.157346935	9.812177306	11.11634396	
ENSG00000001167.13	11.26795708	10.67683861	10.71338651	11.86534664	
ENSG00000001460.16	9.853309555	7.794415866	8.550746785	9.519636253	
ENSG00000001461.15	10.69348696	8.864186145	9.618385502	9.958552715	
ENSG00000001497.15	11.38046107	10.97441459	11.62981194	10.97226185	
```

## 2. Survival Data 
The survival data structure consists of the following format:

First column: Sample IDs (same IDs as the expression matrix columns)
Second column: Indicates survival status (0: alive, 1: dead)
Third column: Additional patient information (e.g., patient ID)
Fourth column: Overall Survival (OS) time in months.
Example:

```
sample			OS	_PATIENT	OS.time
TCGA-VD-A8KM-01A	0	TCGA-VD-A8KM	4
TCGA-VD-AA8R-01A	0	TCGA-VD-AA8R	6
TCGA-VD-AA8M-01A	0	TCGA-VD-AA8M	6
TCGA-WC-A884-01A	0	TCGA-WC-A884	12
TCGA-VD-A8KG-01A	0	TCGA-VD-A8KG	19
TCGA-VD-A8KF-01A	1	TCGA-VD-A8KF	40
TCGA-VD-AA8N-01A	0	TCGA-VD-AA8N	44
TCGA-VD-AA8T-01A	0	TCGA-VD-AA8T	49
TCGA-VD-A8KK-01A	0	TCGA-VD-A8KK	64
TCGA-VD-A8KN-01A	0	TCGA-VD-A8KN	68
TCGA-WC-A88A-01A	1	TCGA-WC-A88A	82
TCGA-VD-AA8P-01A	0	TCGA-VD-AA8P	86
TCGA-VD-A8KD-01A	1	TCGA-VD-A8KD	114
```
## 3. Outputs

 a) partial likelihood deviation found in the 5 fold CV

  ![grafico1LASSOcolorrectalredid](https://github.com/almejiaga/Lassocox/assets/124840761/3ad54acf-55f3-4bfb-a667-d9310d8b64c9)

  b) coefficient pathway for each significant (non-zero gene)

  ![grafico2LASSOcolorrectal](https://github.com/almejiaga/Lassocox/assets/124840761/84fad7ca-2fd2-4031-abfd-7851683bd509)


3. Running the script




