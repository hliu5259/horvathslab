# Seurat1-2.R
# LAST UPDATED ON DEC.22 2020
# filter the Seurat datset based on the feature distribution(saved in the '_beforefilter.rds')

# load package
print('loading required package data.table, tidyverse, Seurat, SingleR')
suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(tidyverse, quietly = TRUE))
suppressMessages(library(Seurat, quietly = TRUE))
suppressMessages(library(SingleR, quietly = TRUE))

#setup graph environment for mgpc
options(bitmapType="cairo")

# manages command line arguments
handle_command_args <- function(args) {
  # make sure all flags are paired
  if(length(args) %%2 != 0) stop("Command line arguments supplied incorrectly!")
  
  # load the flags into a "dictionary"
  arg_df <- data.frame(cbind(flag = args[seq(1, length(args), by = 2)], value = args[seq(2, length(args), by = 2)])) %>% mutate_all(as.character)
    
  # imput sampel list and min feature, max feature, number of mito gene contians
  sample <<- arg_df$value[arg_df$flag == "-s"]
  sample_list <<- fread(sample, header = F)
  sample_id <<- sample_list$V1
  nfeature_min <- as.numeric(sample_list$V2)
  nfeature_max <- as.numeric(sample_list$V3)
  n_percent <- as.numeric(sample_list$V4)
  
  # check the argument
  if (nrow(sample_list) ==0) stop("sample id supplied incorrectly")
  if (ncol(sample_list) < 4) stop("feature argment applied incorrectly")
  
  # check the import dataset
  for (i in 1:nrow(sample_list)){
    if (file.exists(paste0(sample_id[i], "_beforefilter.rds") == F)
            print(paste0(sample_id[i], "_beforefilter.rds NOT exist"))
    else print(paste0("gonna process ", sample_id[i]))
  }
}

# read arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
handle_command_args(args)


for (i in 1:nrow(sample_list)){
    # read Seurat datset
    Gene_Seurat <- readRDS(paste0(sample_id[i],"_beforefilter.rds"))
    # write outputs
    if (!dir.exists(sample_id[i])) {
      cat('Creating output directory...\n')
      dir.create(sample_id[i], showWarnings = FALSE)
      
    # filter unique feature
    Gene_Seurat <- subset(x = Gene_Seurat, subset = nFeature_RNA >nfeature_min[i] & nFeature_RNA <nfeature_max[i] & percent.mt <n_percent[i] )
    
    # plot filtered feature distribution
    png(paste0(sample_id[i],"_feature_distribution_filtered.png"), width = 850, height = 400)
    plot1 <- FeatureScatter(object = Gene_Seurat, feature1 = "nCount_RNA", feature2 = "percent.mt") +
      geom_point()+ scale_fill_viridis_c() + stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
      geom_hex(bins = 70)  + theme_bw()
    plot2 <- FeatureScatter(object = Gene_Seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
      geom_point() + stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +scale_fill_viridis_c() +
      geom_hex(bins = 70)  + theme_bw()
    print(plot1 + plot2)
    dev.off()
    
    # save filterd Seurat rds file
    saveRDS(Gene_Seurat,paste0(sample_id[i],"_filtered.rds"))
    
    # remove the orgrial rds file to release more space
    file.remove(paste0(sample_id[i]), "_beforefilter.rds")
    
    # normatlizaion and scale
    Gene_Seurat <- SCTransform(object = Gene_Seurat, vars.to.regress = "percent.mt", verbose = FALSE, variable.features.n = 6000)
    nor_gene_matrix <- as.data.frame(Gene_Seurat@assays[["RNA"]]@data)
    
    # build up the normatlized gene expression mattrix
    fwrite(nor_gene_matrix, paste0(sample_id[i],"_GE_matrix_filtered.txt"), sep = '\t', quote = F, row.names = T)

    # PCA and cluster
    Gene_Seurat <- RunPCA(Gene_Seurat, verbose = FALSE, npcs = 20)
    
    # plot pca distribution
    png(paste0(sample_id[i],"_Seurat_pca.png"), width = 450, height = 400)
    print(ElbowPlot(object = Gene_Seurat))
    dev.off()
    
    # use UMAP as visualization methods
    Gene_Seurat<- RunUMAP(Gene_Seurat, dims = 1:20)
    Gene_Seurat <- FindNeighbors(Gene_Seurat, verbose = FALSE, dims = 1:20)
    Gene_Seurat <- FindClusters(Gene_Seurat, verbose = FALSE, resolution = 0.2)

    png(paste0(sample_id[i],'_clusters_Seurat_umap.png'), width = 450, height = 400)
    print(DimPlot(Gene_Seurat, label = TRUE, reduction = "umap", group.by ="seurat_clusters") + NoLegend())
    dev.off()
    
    # build up Seurat cluster
    ide <- data.frame(Gene_Seurat@active.ident)
    fwrite(ide, paste0(sample_id[i],'_clusters_Seurat.txt'),sep = "\t",quote = F, row.names = T, col.names = T)

    # use SingleR to annotate cell type
    # set up blueprintEncodeData as reference datset
    rna_re <- BlueprintEncodeData()
    b <- GetAssayData(Gene_Seurat)
    
    # get cluster information from Seurat
    cluster <- Gene_Seurat@active.ident
    
    # link Seurat cluster to SinlgR cluster
    result_cluster <- SingleR(test = b, ref = rna_re, labels = rna_re$label.fine, method="cluster", clusters = cluster)
    Gene_Seurat[["SingleR.cluster.labels"]] <-
      result_cluster$labels[match(Gene_Seurat[[]]["seurat_clusters"]$seurat_clusters, rownames(result_cluster))]
    png(paste0(sample_id[i],"_beforebatchcc_SingleR.png"), width = 450, height = 300)
    print(DimPlot(Gene_Seurat, group.by =  "SingleR.cluster.labels", reduction = "umap", label = TRUE) + labs(title = sample_id[i]))
    dev.off()
    # plot heatmap
    png(paste0(sample_id[i],"_filtered_heatmap.png"), width = 450, height = 300)
    print(plotScoreHeatmap(result_cluster))
    dev.off()
    
    # save Seurat clustered sinleR annotated rds file
    saveRDS(Gene_Seurat, file = paste0(sample_id[i],"_Seurat_clustered_singleR.rds"))


}
#import featureCount gene_expression matrix
