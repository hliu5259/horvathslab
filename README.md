# scReQTL: an approach to correlate SNVs to gene expression from individual scRNA-seq datasets
Using for scReQTL: an approach to correlate SNVs to gene expression from individual scRNA-seq datasets

update by Nov.22 2020

Gene expression matrix analysis for N8,N7,N5 sample from Liu, X., Xiang, Q., Xu, F. et al. Single-cell RNA-seq of cultured human adipose-derived mesenchymal stem cells. Sci Data 6, 190031 (2019). https://doi.org/10.1038/sdata.2019.31

# Seurat Pipeline for analysis GE matrix before intergate
Updated Nov.22 2020

## [Seurat1-1.R](https://github.com/hliu5259/scReQTL/blob/master/seurat1-1.R)
### Command-line:
	Rscript Seurat1-1.R -s <sample_list>

### Description

This script is to generate the Seurat original datset to visualization the feature distribution for downstream analysis. 


### Input 
A list containing the sample_id 
A directory containing the gene expression files (one per sample) 

### sample input matrix file name
	<sample_id>_wide_counts.tsv 

### Required Argument
	-s Sample list contains sample_id
	

## [Seurat1-2.R](https://github.com/hliu5259/scReQTL/blob/master/seurat1-2.R)
### Command-line:
	Rscript Seurat1-2.R -s <sample_feauture_list>

### Description

This script is to filter the Seurat original datset to do cluster.


### Input 
A list containing the sample_id, feature_min, feature_max, percent of Mtio
A directory containing the original Seurat files (_beforefilter.rds)

### sample input matrix file name
	<sample_id>_beforefilter.rds

### Required Argument
	-s Sample_feature_list contains sample_id, feature_min, feature_max, percent of Mito
