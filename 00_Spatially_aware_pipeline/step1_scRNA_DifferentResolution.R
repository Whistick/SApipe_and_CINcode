library(harmony)
library(SCP)
library(Seurat)
library(dplyr)
library(qs)

adata <- qread("adata.qs")
# set different resolutions 
adata <- FindClusters(adata, verbose = FALSE,resolution =0.1)
CellDimPlot(srt = adata, group.by = "seurat_clusters",reduction = "umap",label_point_color="white",label_repel=T,label_segment_color="white", theme_use = "theme_blank",label=T,label_insitu=T,xlab="UMAP-1",ylab="UMAP-2")

# set different resolutions 
adata <- FindClusters(adata, verbose = FALSE,resolution =0.2)
CellDimPlot(srt = adata, group.by = "seurat_clusters",reduction = "umap",label_point_color="white",label_repel=T,label_segment_color="white", theme_use = "theme_blank",label=T,label_insitu=T,xlab="UMAP-1",ylab="UMAP-2")

# set different resolutions 
adata <- FindClusters(adata, verbose = FALSE,resolution =0.3)
CellDimPlot(srt = adata, group.by = "seurat_clusters",reduction = "umap",label_point_color="white",label_repel=T,label_segment_color="white", theme_use = "theme_blank",label=T,label_insitu=T,xlab="UMAP-1",ylab="UMAP-2")


##### Preprocess 
# adata <- NormalizeData(adata, normalization.method = "LogNormalize")
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
adata <- ScaleData(adata,features = rownames(adata))
adata <- RunPCA(adata, features = VariableFeatures(adata), assay = "RNA")
adata <- FindNeighbors(adata, dims = 1:15)
adata <- FindClusters(adata, verbose = FALSE,resolution =0.2)
adata <- RunUMAP(adata, dims = 1:15)
# batch effect 
library(harmony)
library(SCP)
adata <- RunHarmony(adata, "Cohort")
adata <- FindNeighbors(adata, reduction = "harmony", dims = 1:15)
adata <- FindClusters(adata, verbose = FALSE,resolution =0.6)
adata <- RunUMAP(adata, reduction = "harmony",dims=1:15)
CellDimPlot(srt = adata, group.by = "seurat_clusters",reduction = "umap",label_point_color="white",label_repel=T,label_segment_color="white", theme_use = "theme_blank",label=T,label_insitu=T,xlab="UMAP-1",ylab="UMAP-2")
qsave("adata.qs")

# Convert Seurat to H5AD

