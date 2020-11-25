# scReQTL: an approach to correlate SNVs to gene expression from individual scRNA-seq datasets
Using for scReQTL: an approach to correlate SNVs to gene expression from individual scRNA-seq datasets

### R-packages required (dependencies)
	* data.table
	* tidyverse
	* Seurat(v3)
	* biomaRt
	* SingleR

# Sample Quality Assessment and Filtering

## [Seurat1-1.R](https://github.com/hliu5259/scReQTL/blob/master/seurat1-1.R)
### Command-line:
	Rscript Seurat1-1.R -s <sample_list>

### Description

This script visualizes the samples facilitating filtering of unwanted cells. It creates graphs depicting feature distribution for each of the samples. 


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
