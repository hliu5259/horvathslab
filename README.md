# scReQTL: an approach to correlate SNVs to gene expression from individual scRNA-seq datasets
Using for scReQTL: an approach to correlate SNVs to gene expression from individual scRNA-seq datasets

### R-packages required (dependencies)
	* data.table
	* tidyverse
	* Seurat(v3)
	* biomaRt
	* SingleR

# Sample Quality Assessment and Filtering

## [seurat1-1.R](https://github.com/hliu5259/scReQTL/blob/master/seurat1-1.R)
### Command-line:
	Rscript seurat1-1.R -s <sample_list> -p <pattern>

### Description

This script visualizes the samples, facilitating filtering of unwanted cells. It creates graphs depicting features' distribution for each of the samples. 

### Input 

A list containing the sample_id (inputs are assumed to be in the current working directory as).

Example of input file:
`head -3 sample_list.txt
 SRR1592398
 SRR1592399
 SRR1592400` 
 when you run the script as `Rscript seurat1-1.R -s sample_list` the script looks for SRR1592398_wide_counts.tsv, SRR1592398_wide_counts.tsv, etc. if the suffix is different then provide it with `-p <pattern>` optional argument.

### Output

This script produces 2 png files and 1 rds file per sample.
* <sample_name>\_feature\_distribution\_vlnplot.png
* <sample_name>\_feature\_distribution.png
* <sample_name>\_beforefilter.rds
The rds file(s) produced here is needed for the next-script.


## [seurat1-2.R](https://github.com/hliu5259/scReQTL/blob/master/seurat1-2.R)
### Command-line:
	Rscript seurat1-2.R -s <sample_feauture_list>

### Description

For each of the samples processed through seurat1-1.R script this script creates Seurat clusters, assigns cell-types to those clusters through SingleR


### Input 
A list containing the sample_id, feature_min, feature_max, percent of Mtio
A directory containing the original Seurat files (_beforefilter.rds_)

### Output

This script produces 5 png files and 1 rds file per sample.
* <sample_name>\_feature\_distribution\_filtered.png
* <sample_name>\_Seurat\_pca.png
* <sample_name>\_clusters\_Seurat\_umap.png
* <sample_name>\_beforebatchcc\_SingleR.png
* <sample_name>\_filtered\_heatmap.png
* <sample_name>\_Seurat\_clustered\_singleR.rds
The rds file(s) produced here is needed for the next-script.


## [Seurat2_mt_batch_cc.R](https://github.com/hliu5259/scReQTL/blob/master/Seurat2_mt_batch_cc.R)
Note: This script (unlike seurat1-1.R and seurat1-2.R) is hardcoded for N5, N7 and N8 samples (specific to paper) with respect to input and quality-specific metrics (mitochondrial content, features, etc.)

### Description
For the N5, N7, and N8 samples (after being processed through)
