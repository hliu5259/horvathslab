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

This script visualizes the features' distribution for each of the samples, to be used for filtering of unwanted cells. 
### Input 
-s <sample_list> 
: A list containing the sample_id (inputs are assumed to be in the current working directory).

Example:
```
head -3 sample_list.txt
 SRR1592398
 SRR1592399
 SRR1592400
 ```

 When run as `Rscript seurat1-1.R -s sample_list` the script looks for SRR1592398_wide_counts.tsv, SRR1592398_wide_counts.tsv, etc, which are gene-cell matrices with raw counts. Gene-cell matrices in our example are generated using featureCounts. For gene-cell matrices with different suffix please supply it with the `-p <pattern>` optional argument.

### Output

For each sample in the `sample_list` the script produces 2 .png files and 1 .rds file, for example:
* SRR1592398\_feature\_distribution\_vlnplot.png
* SRR1592398\_feature\_distribution.png
* SRR1592398\_beforefilter.rds

The .rds file(s) produced here is needed for the next-script.

## [seurat1-2.R](https://github.com/hliu5259/scReQTL/blob/master/seurat1-2.R)
### Command-line:
	Rscript seurat1-2.R -s <sample_feauture_list>

### Description

For each of the samples processed through seurat1-1.R script seurat1-2.R filters the unwanted cells as specified in the <sample_feauture_list>, and clusters the cells using the SCTransform function of Seurat (v.3.0 or later), after which assigns cell-types to those clusters using SingleR.


### Input 
-s <sample_feauture_list>
: A list containing the sample_id, minimum # of genes/cell (nfeature_min), maximum # of genes/cell (nfeature_max), maximum percentage of mitochondrial genes/cell

Example:
```
 SRR1592398 3000 8000 6
 SRR1592399 2500 7500 5
 SRR1592400 3000 8000 7
```

### Output

This script produces 5 .png files and 1 .rds file per sample.

Example:
* SRR1592398\_feature\_distribution\_filtered.png
* SRR1592398\_Seurat\_pca.png
* SRR1592398\_clusters\_Seurat\_umap.png
* SRR1592398\_beforebatchcc\_SingleR.png
* SRR1592398\_filtered\_heatmap.png
* SRR1592398\_Seurat\_clustered\_singleR.rds

The .rds file(s) produced here is needed for the next-script.

## [Seurat2_mt_batch_cc.R](https://github.com/hliu5259/scReQTL/blob/master/Seurat2_mt_batch_cc.R)
Note: This script (unlike seurat1-1.R and seurat1-2.R) is hardcoded for N5, N7 and N8 samples (specific to our example analysis) with respect to input and quality-specific metrics (mitochondrial content, features, etc.)

### Description
The N5, N7, and N8 samples, after being processed through seurat1-1.R and seurat1-2.R (the `rds` files), are used as inputs here. Mitochondrial genes are filtered out completely.
This script performs 2 main jobs:

1. Integrates the samples (using SCTransform function), removes technical variations, plots the feature distributions, plots the UMAP clusters, and annotates and plots cell-type information through SingleR.
2. The integrated dataset is divided into individual samples and for each of them, cell-cycle effects are regressed out, PCs are analyzed (JackStraw and Elbow plots), cell-types are annotated through SingleR and plotted; the cluster information is stored. 

### Output

This script produces 4 .png files (integrated for all the samples), 7 .png files per sample, 2 .txt files per sample, and 1 .rds file per sample.

**Integrated Data**
* N578\_feature\_distribution\_vlnplot.png
* N578\_clusters\_Seurat\_umap.png
* N578\_clusters\_SingleR\_umap.png
* N578\_clusters\_SingleR\_sample\_split\_umap.png


**Per-sample**
* <sample_name>\_feature\_distribution\_filtered.png
* <sample_name>\_Seurat\_pca.png
* <sample_name>\_clusters\_Seurat\_umap.png
* <sample_name>\_beforebatchcc\_SingleR.png
* <sample_name>\_filtered\_heatmap.png
* <sample_name>\_Seurat\_clustered\_singleR.rds

The rds file(s) produced here is needed for the next-script.
