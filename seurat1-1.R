# Seurat1-1.R
# LAST UPDATED ON DEC.22 2020
# This script is to generate the Seurat original datset to visualization the feature distribution for downstream analysis

# load package
suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(tidyverse, quietly = TRUE))
suppressMessages(library(Seurat, quietly = TRUE))
suppressMessages(library(biomaRt, quietly = TRUE))

# manages command line arguments
handle_command_args <- function(args) {
  # make sure all flags are paired
  if(length(args) %%2 != 0) stop("Command line arguments supplied incorrectly!")
  
  # load the flags into a "dictionary"
  arg_df <- data.frame(cbind(flag = args[seq(1, length(args), by = 2)], value = args[seq(2, length(args), by = 2)])) %>% mutate_all(as.character)
    
  # sample list
  sample <<- arg_df$value[arg_df$flag == "-s"]
  sample_list <<- fread(sample, header = F)
  sample_id <<- sample_list$V1
  
  if (nrow(sample_list) ==0) stop('sample id supplied incorrectly')
 
  # check the import dataset
  for (i in 1:nrow(sample_list)){
    if (file.exists(paste0(sample_id[i], '_wide_counts.tsv') == F) print(paste0('provided ', sample_id[i],
                                                        ' NOT exist'))
    else { print(paste0('gonna process ', sample_id[i]))}
  }
}
# read arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
handle_command_args(args)

for (i in 1:nrow(sample_list)){
    #import featureCount gene_expression matrix
    Gene_matrix <- fread(paste0(sample_id[i], '_wide_counts.tsv'))
    #Gene_matrix <- spread(Gene_matrix, cell, count)
    # ensemble id change to gene_name
    library(biomaRt)
    ##listEnsemblArchives()
    # use human sapines ensembl verson100 data to annotate (should change by different dataset, keep the lastest version)
    ensemble <- useEnsembl(dataset = "hsapiens_gene_ensembl", biomart = "ensembl", version =100)
    query <- getBM(attributes =
                     c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position"),
                   filters = "ensembl_gene_id", values = Gene_matrix$gene, mart = ensemble)
    colnames(query)[1] <- "gene"
    Gene_matrix_name <- inner_join(query,Gene_matrix, by = "gene")
    Gene_matrix_name <- Gene_matrix_name[,-c(1,3,4,5)]
    # deduplicate file
    index <- duplicated(Gene_matrix_name$external_gene_name)
    Gene_matrix_name <- Gene_matrix_name[!index,]
    rownames(Gene_matrix_name) <- Gene_matrix_name$external_gene_name
    Gene_matrix_name <- Gene_matrix_name[,-1]
    
    #import dataset into Seurat as min.cells = 200, min.features =2000
    Gene_Seurat <- CreateSeuratObject(counts = Gene_matrix_name ,project = sample_id , min.cells = 300, min.features =2000)
    #calculate mito QC metrics
    Gene_Seurat <- PercentageFeatureSet(Gene_Seurat, pattern = "^MT-", col.name = "percent.mt")
    # visualization for the feature distribution
    png(paste0(sample_id,'_','feature_distribution_vlnplot.png'), width = 850, height = 400)
    print(VlnPlot(object = Gene_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
    dev.off()
    # density plot of feature feature relationship
    png(paste0(sample_id, '_feature_distribution.png'), width = 850, height = 400)
    plot1 <- FeatureScatter(object = Gene_Seurat, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_point()+ scale_fill_viridis_c() + stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") + geom_hex(bins = 70)  + theme_bw()
    plot2 <- FeatureScatter(object = Gene_Seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ geom_point() + stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +scale_fill_viridis_c() + geom_hex(bins = 70)  + theme_bw()
    print(plot1 + plot2)
    dev.off()
    # saving original Seurat dataset
    saveRDS(Gene_Seurat,paste0(sample_id[i],'_beforefilter.rds'))
}

print('DONE')
