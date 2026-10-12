library(spatstat.explore)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

rm(list=ls())

options(future.globals.maxSize= 891289600000)

SAMPLENAME=c("J259B","J259X")

data=list() # data list


### read in the data
for(sInd in 1:length(SAMPLENAME)){
  print(c(sInd, SAMPLENAME[sInd]))
  ## read in data
  data_dir <- paste0("../code/",SAMPLENAME[sInd],"/outs/")
  list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
  data[[sInd]] <- Load10X_Spatial(data.dir = data_dir)
  
  data[[sInd]]$library=SAMPLENAME[sInd]
  
  print(data[[sInd]])
}
names(data)=SAMPLENAME

# [1] "1"     "J259B"
# An object of class Seurat 
# 18085 features across 2397 samples within 1 assay 
# Active assay: Spatial (18085 features, 0 variable features)
# 1 image present: slice1
# [1] "2"     "J259X"
# An object of class Seurat 
# 18085 features across 3923 samples within 1 assay 
# Active assay: Spatial (18085 features, 0 variable features)
# 1 image present: slice1


####### Data pre-processing: 
for(sInd in 1:length(data)){
  print(c(sInd, SAMPLENAME[sInd]))
  
  ## data visualization
  pdf(paste0("Vln.",sInd,".",SAMPLENAME[sInd],".pdf"), width=10, height=5)
  plot1 <- VlnPlot(data[[sInd]], features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  plot2 <- SpatialFeaturePlot(data[[sInd]], features = "nCount_Spatial") + theme(legend.position = "right")
  print(wrap_plots(plot1, plot2))
  dev.off()
  
  ## filter out zero cells (add this step, incase SCTransform will show error)
  # data[[i]] <- subset(x = data[[i]], subset = nCount_Spatial > 0)
  # print(data[[i]])
  
  ## transformation
  data[[sInd]] <- SCTransform(data[[sInd]], assay = "Spatial", verbose = FALSE)
  
  ## Gene expression visualization
  # ALLGENE=rownames(data[[i]]@assays$SCT)
  # plotGene=c("Alb","Gapdh")
  # plotGene%in%ALLGENE
  # 
  # pdf(paste0("Feature.",i,".pdf"), width=10, height=5)
  # print(SpatialFeaturePlot(data[[i]], features = plotGene))
  # dev.off()
  
}

save(data, file="data.rdata")

load("data.rdata")


#### integration, pipeline 1, as suggested by ST pipeline, but no normalization
# OBJ.merge=merge(x=data[[1]], y=data[[2]], ## list all
#                 add.cell.ids=SAMPLENAME,
#                 project="Spatial")
# 
# OBJ.merge
# # An object of class Seurat 
# # 36126 features across 6320 samples within 2 assays 
# # Active assay: SCT (18041 features, 0 variable features)
# # 1 other assay present: Spatial
# # 2 images present: slice1, slice1_J259X
# 
# DefaultAssay(OBJ.merge) <- "SCT"
# 
# tmp=c()
# for(i in 1:length(data)){
#   tmp=c(tmp, VariableFeatures(data[[i]]))
# }
# VariableFeatures(OBJ.merge) <- tmp
# 
# OBJ.merge <- RunPCA(OBJ.merge, verbose = FALSE)
# OBJ.merge <- FindNeighbors(OBJ.merge, dims = 1:30)
# OBJ.merge <- FindClusters(OBJ.merge, verbose = FALSE, resolution=0.8)
# table(OBJ.merge@meta.data$seurat_clusters)
# 
# OBJ.merge <- RunUMAP(OBJ.merge, dims = 1:30)
# 
# ## add in the library information
# LIB=sapply(strsplit(rownames(OBJ.merge@meta.data), split="_"), function(x) return(x[1]))
# OBJ.merge@meta.data$library=LIB
# 
# table(OBJ.merge@meta.data$library)
# # J259B J259X 
# # 2397  3923 
# 
# saveRDS(file = "OBJ.merge.RDS", OBJ.merge)
# 
# OBJ.merge <- readRDS("OBJ.merge.RDS")
# 
# 
# table(OBJ.merge@meta.data$seurat_clusters, OBJ.merge@meta.data$library )
# 
# # J259B J259X
# # 0      0  1107
# # 1      9   899
# # 2      6   822
# # 3    685     2
# # 4    656    14
# # 5    385    85
# # 6    405     0
# # 7      5   378
# # 8      0   325
# # 9    246     0
# # 10     0   150
# # 11     0   141




############# pipeline 2: SCTransform, as scRNA-seq data

### Preform the Data integration between two samples using 3000 features

select.features <- SelectIntegrationFeatures(object.list = data, nfeatures = 3000)

data <- PrepSCTIntegration(object.list = data,
                                anchor.features = select.features, 
                                verbose = TRUE)


select.anchors <- FindIntegrationAnchors(object.list = data, 
                                       normalization.method = "SCT", 
                                       anchor.features = select.features,
                                       verbose = T,reference = 1,)

OBJ.merge <- IntegrateData(anchorset = select.anchors, normalization.method = "SCT", 
                                 verbose = FALSE,)


## Generate the PCA, UMAP and visualizations
OBJ.merge <- RunPCA(OBJ.merge, verbose = FALSE)
OBJ.merge <- RunUMAP(OBJ.merge, dims = 1:20,n.neighbors = 30)

DimPlot(OBJ.merge, label=T,group.by = c("library"), combine = FALSE)


##Generate initial cluster calls
OBJ.merge <- FindNeighbors(OBJ.merge, reduction = "pca", dims = 1:20)
OBJ.merge <- FindClusters(OBJ.merge, resolution = 1) ## old 9 clusters with resolution=0.5

### swap and scale to RNA assay for visuals
DefaultAssay(OBJ.merge) <- "Spatial"
OBJ.merge <- FindVariableFeatures(OBJ.merge)
OBJ.merge <- NormalizeData(OBJ.merge)
OBJ.merge <- ScaleData(OBJ.merge)

FeaturePlot(OBJ.merge,features=c("ALB","GAPDH"),split.by="library")

names(OBJ.merge@images)=SAMPLENAME

saveRDS(OBJ.merge, file = "OBJ.merge.RDS")
# OBJ.merge=readRDS("OBJ.merge.RDS")

### plot

pdf(paste0("DimPlot.pdf"), width=10, height=5)
DimPlot(OBJ.merge, reduction = "umap", group.by = c("seurat_clusters", "library"))
dev.off()

pdf(paste0("DimPlot.cluster.pdf"), width=6, height=5)
DimPlot(OBJ.merge, reduction = "umap", group.by = c("seurat_clusters"))
dev.off()

pdf(paste0("SpatialDimPlot.pdf"), width=10, height=5)
SpatialDimPlot(OBJ.merge)
dev.off()

pdf(paste0("SpatialDimPlotUpdate.pdf"), width=10, height=5)
SpatialDimPlot(OBJ.merge, pt.size.factor=80, image.alpha=0.3)
dev.off()


####### check gene expression
DefaultAssay(OBJ.merge) <- "Spatial"

ALLGENE=rownames(OBJ.merge@assays$SCT)

plotGene=c("GLUL","ALB") 
plotGene[!plotGene%in%ALLGENE]

pdf(paste0("SpatialFeaturePlot.pdf"), width=15, height=15)
SpatialFeaturePlot(OBJ.merge, features = plotGene)
dev.off()


table(OBJ.merge@meta.data$seurat_clusters, OBJ.merge@meta.data$library)

# J259B J259X
# 0     16   718
# 1     15   652
# 2     17   627
# 3    209   372
# 4    544     4
# 5     32   429
# 6    350    62
# 7    155   254
# 8    363    26
# 9    323    54
# 10    14   299
# 11    14   198
# 12   197     8
# 13    20   128
# 14    67    73
# 15    61    19



########## marker per cluster
sobj.markers <- FindAllMarkers(OBJ.merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

marker.output.top <- sobj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
print(marker.output.top)

marker.output.all <- sobj.markers %>% group_by(cluster) 
dim(marker.output.all) # 16593     7
table(marker.output.all$cluster)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 328  632  160  188  900  591  415  218 5522  713 1082  449  787  321 2274 2013 
write.csv(marker.output.all, file="marker_for_cluster.csv",row.names = F, quote=F)




####### check gene expression, spatial

DefaultAssay(OBJ.merge) <- "Spatial"

ALLGENE=rownames(OBJ.merge@assays$SCT)


# Hepatocytes
# ALB, TF, TTR, CYP3A4, APOA1
plotGene=c("ALB", "TF", "TTR", "CYP3A4", "APOA1") 
plotGene[!plotGene%in%ALLGENE]

pdf(paste0("SpatialFeaturePlot.Hepatocytes.pdf"), width=10, height=20)
SpatialFeaturePlot(OBJ.merge, features = plotGene, pt.size.factor=50)
dev.off()

pdf(paste0("FeaturePlot.Hepatocytes.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene, ncol=2)
dev.off()


# Liver sinusoidal ECs
# LYVE1, ENG, STAB1, PECAM1, CD36
plotGene=c("LYVE1", "ENG", "STAB1", "PECAM1", "CD36") 
plotGene[!plotGene%in%ALLGENE]

pdf(paste0("SpatialFeaturePlot.SinusoidalEC.pdf"), width=10, height=20)
SpatialFeaturePlot(OBJ.merge, features = plotGene, pt.size.factor=50)
dev.off()

pdf(paste0("FeaturePlot.SinusoidalEC.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene)
dev.off()


# Vascular ECs
# PECAM1, CCL21, MMRN1, GNG11, FLT4
plotGene=c("PECAM1", "CCL21", "MMRN1", "GNG11", "FLT4") 
plotGene[!plotGene%in%ALLGENE]

pdf(paste0("SpatialFeaturePlot.VascularEC.pdf"), width=10, height=20)
SpatialFeaturePlot(OBJ.merge, features = plotGene, pt.size.factor=50)
dev.off()

pdf(paste0("FeaturePlot.VascularEC.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene)
dev.off()


# Cholangiocytes
# EPCAM, KRT7, KRT19, CLDN4, TACSTD2
plotGene=c("EPCAM", "KRT7", "KRT19", "CLDN4", "TACSTD2") 
plotGene[!plotGene%in%ALLGENE]

pdf(paste0("SpatialFeaturePlot.Cholangiocytes.pdf"), width=10, height=20)
SpatialFeaturePlot(OBJ.merge, features = plotGene, pt.size.factor=50)
dev.off()

pdf(paste0("FeaturePlot.Cholangiocytes.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene)
dev.off()


# Stellate cells
# ACTA2, DCN, BGN, COLEC11, IGFBP3
plotGene=c("ACTA2", "DCN", "BGN", "COLEC11", "IGFBP3") 
plotGene[!plotGene%in%ALLGENE]

pdf(paste0("SpatialFeaturePlot.StellateCells.pdf"), width=10, height=20)
SpatialFeaturePlot(OBJ.merge, features = plotGene, pt.size.factor=50)
dev.off()

pdf(paste0("FeaturePlot.StellateCells.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene)
dev.off()


# Kupffer cells
# CD68, CD163, LYZ, C1QA, AIF1
plotGene=c("CD68", "CD163", "LYZ", "C1QA", "AIF1") 
plotGene[!plotGene%in%ALLGENE]

pdf(paste0("SpatialFeaturePlot.KupfferCells.pdf"), width=10, height=20)
SpatialFeaturePlot(OBJ.merge, features = plotGene, pt.size.factor=50)
dev.off()

pdf(paste0("FeaturePlot.KupfferCells.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene)
dev.off()



# T cells
# CD3D, CD2, IL7R, TRBC2, CD69
plotGene=c("CD3D", "CD2", "IL7R", "TRBC2", "CD69") 
plotGene[!plotGene%in%ALLGENE]

pdf(paste0("SpatialFeaturePlot.Tcells.pdf"), width=10, height=20)
SpatialFeaturePlot(OBJ.merge, features = plotGene, pt.size.factor=50)
dev.off()

pdf(paste0("FeaturePlot.Tcells.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene)
dev.off()


# B cells / Plasma cells
# IGKC, JCHAIN, CD79A, CD27, CD74
plotGene=c("IGKC", "JCHAIN", "CD79A", "CD27", "CD74") 
plotGene[!plotGene%in%ALLGENE]

pdf(paste0("SpatialFeaturePlot.Bcells.pdf"), width=10, height=20)
SpatialFeaturePlot(OBJ.merge, features = plotGene, pt.size.factor=50)
dev.off()

pdf(paste0("FeaturePlot.Bcells.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene)
dev.off()


# CD genes and NK cell markers
# Check CD4, CD8, CD11
# NK cell markers: Nkg7, Klrd1, Prf1, Cd7, Trdc

## CD11 -> ITGAM

plotGene=c("CD4","CD8A","ITGAM",
           "NKG7", "KLRD1", "PRF1", "CD7", "TRDC") 
plotGene[!plotGene%in%ALLGENE]

pdf(paste0("SpatialFeaturePlot.CDgene.NKcell.pdf"), width=10, height=20)
SpatialFeaturePlot(OBJ.merge, features = plotGene, pt.size.factor=50)
dev.off()

pdf(paste0("FeaturePlot.CDgene.NKcell.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene)
dev.off()



########### dot plot

plotGene=c("ALB", "TF", "TTR", "CYP3A4", "APOA1",
           "LYVE1", "ENG", "STAB1", "PECAM1", "CD36",
           "PECAM1", "CCL21", "MMRN1", "GNG11", "FLT4",
           "EPCAM", "KRT7", "KRT19", "CLDN4", "TACSTD2",
           "ACTA2", "DCN", "BGN", "COLEC11", "IGFBP3",
           "CD68", "CD163", "LYZ", "C1QA", "AIF1",
           "CD3D", "CD2", "IL7R", "TRBC2", "CD69",
           "IGKC", "JCHAIN", "CD79A", "CD27", "CD74")
plotGene=unique(plotGene)
plotGene[!plotGene%in%ALLGENE]

pdf(paste0("DotPlot.geneMarker.pdf"), width=8, height=15)
DotPlot(object=OBJ.merge, features=plotGene) +
  coord_flip()
dev.off()





