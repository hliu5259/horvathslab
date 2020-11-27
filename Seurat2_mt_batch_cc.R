# Seurat_mt_batch_cc.R
# Last updated on NOV.23 2020
# used for scReQTL (N8, N7,N5 dataset)

# load package
suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(tidyverse, quietly = TRUE))
suppressMessages(library(Seurat, quietly = TRUE))
suppressMessages(library(SingleR, quietly = TRUE))

if (packageVersion("Seurat") < "3.0.0") {
  stop(paste0("You have Seurat version", packageVersion("Seurat"), "installed. Please make sure you have Seurat version > 3.0.0"))
}

# Load all data-----------------------------------------------------------------------------------------
#For these samples: 3k<nfeatures<8k and mt_percent<7
# load filtered Seurat data
N8 <- readRDS('N8_Seurat_clustered_singleR.rds')
N7 <- readRDS('N7_Seurat_clustered_singleR.rds')
N5 <- readRDS('N5_Seurat_clustered_singleR.rds')

# cell cycle calculation
s.genes <- fread('stage_s.txt')
s.genes <- s.genes$Gene
g2m.genes <- fread('stage_G2M.txt')
g2m.genes <- g2m.genes$Gene
#-------------------------------------------------------------------------------------------------------

# Data pre-processing-----------------------------------------------------------------------------------
# remove mitochondrial genes
N8_outMT <- N8@assays[['RNA']]@counts
N8_outMT <- N8_outMT[-grep(pattern = '^MT', row.names(N8_outMT)),  ]
N8 <- CreateSeuratObject(counts = N8_outMT,project = "N8")
rm(N8_outMT)
N7_outMT <- N7@assays[['RNA']]@counts
N7_outMT <- N7_outMT[-grep(pattern = '^MT', row.names(N7_outMT)),  ]
N7 <- CreateSeuratObject(counts = N7_outMT,project = "N7")
rm(N7_outMT)
N5_outMT <- N5@assays[['RNA']]@counts
N5_outMT <- N5_outMT[-grep(pattern = '^MT', row.names(N5_outMT)),  ]
N5 <- CreateSeuratObject(counts = N5_outMT,project = "N5")
rm(N5_outMT)

# normalization and scale the dataset (only use nUMI to regress which is default)
# Set 3000 variable features(default: variable.features.n = 3000)
## our data has been regressed (nUMI) based on umi_tools
N8 <- SCTransform(N8)
N7 <- SCTransform(N7)
N5 <- SCTransform(N5)
#-------------------------------------------------------------------------------------------------------

#integrate three samples together
list <- list(N5, N7, N8)
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 8000)
options(future.globals.maxSize = 9768*1024^2)
list <- PrepSCTIntegration(object.list = list, anchor.features = features,
                           verbose = FALSE)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                  anchor.features = features, verbose = FALSE)
N578 <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)

png("N578_feature_distribution_vlnplot.png", width = 850, height = 400)
VlnPlot(object = N578, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

saveRDS(N578, 'N578_outMT_perPCA.rds')

# cluster
DefaultAssay(N578) <- "integrated"
N578 <- RunPCA(N578, verbose = FALSE)
N578 <- RunUMAP(N578, verbose = FALSE, dim = 1:30)
N578<- FindNeighbors(N578, dims = 1:30) #,  k.param = 10
N578<- FindClusters(N578, resolution = 0.2)

# plot cluster information
png("N578_clusters_Seurat_umap.png", width = 850, height = 400)
DimPlot(N578, group.by = "orig.ident",reduction = "umap", label = TRUE)
dev.off()

# annotation cell type
rna_re <- BlueprintEncodeData()
cluster <- N578@active.ident
b <- GetAssayData(N578)
result_cluster <- SingleR(test = b, ref = rna_re, labels = rna_re$label.fine, method="cluster", clusters = cluster)
N578[["SingleR.cluster.labels"]] <-
  result_cluster$labels[match(N578[[]]["seurat_clusters"]$seurat_clusters, rownames(result_cluster))]
png("N578_clusters_Seurat_umap.png", width = 850, height = 400)
DimPlot(N578, group.by =  "SingleR.cluster.labels", reduction = "umap", label = TRUE)
dev.off()
png("N578_clusters_Seurat_umap.png", width = 850, height = 400)
DimPlot(N578, group.by= "SingleR.cluster.labels", split.by = "orig.ident",reduction = "umap", label = TRUE) + NoLegend()
dev.off()

#split the file
list_new <- SplitObject(N578, split.by = "orig.ident")
N8 <- list_new$N8
N7 <- list_new$N7
N5 <- list_new$N5


# Identifying clusters individually after batch effect removal
# N8-----------------------------------------------------------------------------------------------
DefaultAssay(N8) <- "integrated"
N8 <- CellCycleScoring(N8, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) #  (use head(N8[[]]) to check)
N8 <- RunPCA(N8)
N8 <- RunUMAP(N8, dims = 1:30)
png("N578_clusters_Seurat_umap.png", width = 850, height = 400)
DimPlot(N8, group.by = "Phase", reduction = "umap",label = TRUE)
dev.off()
N8 <- ScaleData(N8, vars.to.regress = c("S.Score", "G2M.Score"),
                 use.umi = FALSE, do.scale = FALSE, do.center =  FALSE )
N8 <- RunPCA(N8)
N8<- RunUMAP(N8, dims = 1:30)#, spread = 1, min.dist = 0.0001, n.epochs = 200, n.neighbors = 10)

png("N578_clusters_Seurat_umap.png", width = 850, height = 400)
DimPlot(N8, group.by = "Phase", reduction = "umap",label = TRUE)
dev.off()

N8<- FindNeighbors(N8, dims = 1:30)#, k.param = 10)
N8<- FindClusters(N8, resolution = 0.2)

png("N578_clusters_Seurat_umap.png", width = 850, height = 400)
ElbowPlot(N8)
dev.off()

N8 <- JackStraw(N8, num.replicate = 100)
N8 <- ScoreJackStraw(N8, dims = 1:20)

png("N578_clusters_Seurat_umap.png", width = 850, height = 400)
JackStrawPlot(N8, dims = 1:20)
dev.off()


# annotation
b <- GetAssayData(N8)
res_re <- BlueprintEncodeData()
cluster <- N8@active.ident
result_cluster <- SingleR(test = b, ref = res_re, labels = res_re$label.fine, method="cluster", clusters = cluster)

png("N578_clusters_Seurat_umap.png", width = 850, height = 400)
plotScoreHeatmap(result_cluster)
dev.off()

N8[["SingleR.cluster.labels"]] <-
  result_cluster$labels[match(N8[[]]["seurat_clusters"]$seurat_clusters, rownames(result_cluster))]

png("N578_clusters_Seurat_umap.png", width = 850, height = 400)
DimPlot(N8, group.by =  "SingleR.cluster.labels", reduction = "umap", label = TRUE)
dev.off()

png("N578_clusters_Seurat_umap.png", width = 850, height = 400)
DoHeatmap(subset(N8, downsample = 100), features = features, size = 3, group.by = "SingleR.cluster.labels")
dev.off()

# write cluster information
ide <- as.data.frame(N8@active.ident)
fwrite(ide, "N8_cluster_sample.txt", row.names = T, sep = '\t', quote = F)
cluster <- data.frame(N8[["SingleR.cluster.labels"]])
fwrite(cluster, "N8_cluster_sample_named.txt", row.names = T, sep = ',', quote = F)

# pca used for N8
N8_pca <- t(as.data.frame(N8@reductions[["pca"]]@cell.embeddings))
fwrite(N8_pca, "N8_pca_matrix.txt", row.names = T, sep = '\t', quote = F)
#-------------------------------------------------------------------------------------------------------

#N7----------------------------------------------------------------------------------------------------
DefaultAssay(N7) <- "integrated"
N7 <- CellCycleScoring(N7, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) #  (use head(N7[[]]) to check)
N7 <- RunPCA(N7)
N7 <- RunUMAP(N7, dims = 1:30)
png("N578_clusters_Seurat_umap.png", width = 850, height = 400)
DimPlot(object = N7, group.by = "Phase", reduction = "umap",label = TRUE)
dev.off()

# Sale N7 sample by the score of S and G2M
N7 <- ScaleData(N7, vars.to.regress = c("S.Score", "G2M.Score"),
                use.umi = FALSE, do.scale = FALSE, do.center =  FALSE )
N7 <- RunPCA(N7)
N7<- RunUMAP(N7, dims = 1:30)#, spread = 1, min.dist = 0.0001, n.epochs = 200, n.neighbors = 10)

png("N578_clusters_Seurat_umap.png", width = 850, height = 400)
DimPlot(object = N7, group.by = "Phase", reduction = "umap",label = TRUE)
dev.off()

N7<- FindNeighbors(N7, dims = 1:30)#, k.param = 10)
N7<- FindClusters(N7, resolution = 0.2)

# write the cluster information
ide <- as.data.frame(N7@active.ident)
fwrite(ide, "N7_cluster_sample.txt", row.names = T, sep = '\t', quote = F)
cluster <- data.frame(N7[["SingleR.cluster.labels"]])
fwrite(cluster, "N7_cluster_sample_named.txt", row.names = T, sep = ',', quote = F)

# annotation
b <- GetAssayData(N7)
cluster <- N7@active.ident
result_cluster <- SingleR(test = b, ref = res_re, labels = res_re$label.fine, method="cluster", clusters = cluster)
N7[["SingleR.cluster.labels"]] <-
  result_cluster$labels[match(N7[[]]["seurat_clusters"]$seurat_clusters, rownames(result_cluster))]

png("N578_clusters_Seurat_umap.png", width = 850, height = 400)
DimPlot(N7, group.by =  "SingleR.cluster.labels", reduction = "umap", label = TRUE)
dev.off()

png("N578_clusters_Seurat_umap.png", width = 850, height = 400)
DoHeatmap(subset(N7, downsample = 100), features = features, size = 3, group.by = "SingleR.cluster.labels")
dev.off()

# pca used for N7
N7_pca <- t(as.data.frame(N7@reductions[["pca"]]@cell.embeddings))
fwrite(N7_pca, "N7_pca_matrix.txt", row.names = T, sep = '\t', quote = F)
#-------------------------------------------------------------------------------------------------------

#N5-----------------------------------------------------------------------------------------------------
DefaultAssay(N5) <- "integrated"
N5 <- CellCycleScoring(N5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) # #  (use head(N5[[]]) to check)
N5 <- RunPCA(N5)
N5 <- RunUMAP(N5, dims = 1:30)

png("N5_beforess_umap.png", width = 550, height = 500)
DimPlot(object = N5, group.by = "Phase", reduction = "umap",label = TRUE)
dev.off()

# Scale N5 sample by the score of S and G2M
N5 <- ScaleData(N5, vars.to.regress = c("S.Score", "G2M.Score"),
                use.umi = FALSE, do.scale = FALSE, do.center =  FALSE )
N5 <- RunPCA(N5)
N5<- RunUMAP(N5, dims = 1:30)#, spread = 1, min.dist = 0.0001, n.epochs = 200, n.neighbors = 10)

png("N5_afterss_umap.png", width = 550, height = 500)
DimPlot(object = N5, group.by = "Phase", reduction = "umap",label = TRUE)
dev.off()

N5<- FindNeighbors(N5, dims = 1:30)#, k.param = 10)
N5<- FindClusters(N5, resolution = 0.2)

# write cluster information
ide <- as.data.frame(N5@active.ident)
fwrite(ide, "N5_cluster_sample.txt", row.names = T, sep = '\t', quote = F)
cluster <- data.frame(N5[["SingleR.cluster.labels"]])
fwrite(cluster, "N5_cluster_sample_named.txt", row.names = T, sep = ',', quote = F)

# annotation
b <- GetAssayData(N5)
cluster <- N5@active.ident
result_cluster <- SingleR(test = b, ref = res_re, labels = res_re$label.fine, method="cluster", clusters = cluster)
N5[["SingleR.cluster.labels"]] <-
  result_cluster$labels[match(N5[[]]["seurat_clusters"]$seurat_clusters, rownames(result_cluster))]

png("N578_clusters_Seurat_umap.png", width = 850, height = 400)
DimPlot(N5, group.by =  "SingleR.cluster.labels", reduction = "umap", label = TRUE)
dev.off()

png("N578_clusters_Seurat_umap.png", width = 850, height = 400)
DoHeatmap(subset(N5, downsample = 100), features = features, size = 3, group.by = "SingleR.cluster.labels")
dev.off()

# pca used for N5
N5_pca <- t(as.data.frame(N5@reductions[["pca"]]@cell.embeddings))
fwrite(N5_pca, "N5_pca_matrix.txt", row.names = T, sep = '\t', quote = F)
#-------------------------------------------------------------------------------------------------------

# save output
saveRDS(N8, "N8_ccregress.rds")
saveRDS(N7, "N7_ccregress.rds")
saveRDS(N5, "N5_ccregress.rds")
