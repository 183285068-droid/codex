
## user R version 4.2.0 and Seurat V4 version


library(spatstat.explore)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

rm(list=ls())

options(future.globals.maxSize= 891289600000)

SAMPLENAME=c("J259B","J259X")

OBJ.merge=readRDS("OBJ.merge.RDS")

SpatialDimPlot(OBJ.merge)


### Define immune spots: average expression of T/B/Kupffer/NK cells

ALLGENE=rownames(OBJ.merge@assays$SCT)

# Kupffer cells
# CD68, CD163, LYZ, C1QA, AIF1
plotGene=c("CD68", "CD163", "LYZ", "C1QA", "AIF1") 
plotGene[!plotGene%in%ALLGENE]
Kupffer.exp=apply(as.matrix(GetAssayData(object = OBJ.merge, slot = "data")[plotGene,]), 2, mean)

# T cells
# CD3D, CD2, IL7R, TRBC2, CD69
plotGene=c("CD3D", "CD2", "IL7R", "TRBC2", "CD69") 
plotGene[!plotGene%in%ALLGENE]
Tcell.exp=apply(as.matrix(GetAssayData(object = OBJ.merge, slot = "data")[plotGene,]), 2, mean)

# B cells / Plasma cells
# IGKC, JCHAIN, CD79A, CD27, CD74
plotGene=c("IGKC", "JCHAIN", "CD79A", "CD27", "CD74") 
plotGene[!plotGene%in%ALLGENE]
Bcell.exp=apply(as.matrix(GetAssayData(object = OBJ.merge, slot = "data")[plotGene,]), 2, mean)

# CD genes and NK cell markers
# Check CD4, CD8, CD11
# NK cell markers: Nkg7, Klrd1, Prf1, Cd7, Trdc
## CD11 -> ITGAM
plotGene=c("CD4","CD8A","ITGAM",
           "NKG7", "KLRD1", "PRF1", "CD7", "TRDC") 
plotGene[!plotGene%in%ALLGENE]
NKcell.exp=apply(as.matrix(GetAssayData(object = OBJ.merge, slot = "data")[plotGene,]), 2, mean)


immune.exp=cbind(Kupffer.exp, Tcell.exp, Bcell.exp, NKcell.exp)



###  check correlation
library(corrplot)
library(RColorBrewer)
M <-cor(immune.exp)
M
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))

#### get all high
immune.ave=apply(immune.exp, 1, mean)
sum(immune.ave>=1) # 77
sum(immune.ave>=1)/length(immune.ave) # 0.01218354

sum(immune.ave>=0.9) # 289
sum(immune.ave>=0.9)/length(immune.ave) # 0.04572785

sum(immune.ave>=0.8) # 771
sum(immune.ave>=0.8)/length(immune.ave) # 0.1219937


ImmuneLabel=ifelse(immune.ave>=0.9, "Immune","NonImmune")
table(ImmuneLabel)
# Immune NonImmune 
# 289      6031 


## select the surrounding cells of immune cells

# J259B, tissue=1
posFile=OBJ.merge@images$J259B@coordinates
inCell=intersect(rownames(posFile), names(which(ImmuneLabel=="Immune")))
length(inCell) # 165

surroundCell=c()
for(cellName in inCell){
  surroundCell=c(surroundCell,rownames(posFile)[which( abs(posFile$row-posFile[cellName,"row"])<=1 & abs(posFile$col-posFile[cellName,"col"])<=2 )])
}

surroundCell=unique(surroundCell)
surroundCell=setdiff(surroundCell, inCell)
length(surroundCell) # 213
ImmuneLabel[surroundCell]="SurroundImmune"



# J259X, tissue=2
posFile=OBJ.merge@images$J259X@coordinates
inCell=intersect(rownames(posFile), names(which(ImmuneLabel=="Immune")))
length(inCell) # 124

surroundCell=c()
for(cellName in inCell){
  surroundCell=c(surroundCell,rownames(posFile)[which( abs(posFile$row-posFile[cellName,"row"])<=1 & abs(posFile$col-posFile[cellName,"col"])<=2 )])
}

surroundCell=unique(surroundCell)
surroundCell=setdiff(surroundCell, inCell)
length(surroundCell) # 457
ImmuneLabel[surroundCell]="SurroundImmune"

table(ImmuneLabel)
# Immune      NonImmune SurroundImmune 
# 289           5361            670 

OBJ.merge@meta.data=cbind(OBJ.merge@meta.data, ImmuneLabel=ImmuneLabel)
table(OBJ.merge@meta.data$ImmuneLabel)
# Immune      NonImmune SurroundImmune 
# 289           5361            670 


## plot the figure

pdf("SpatialDimPlot.Immune.Update.pdf", width=12, height=5)
SpatialDimPlot(OBJ.merge, group.by = "ImmuneLabel", pt.size.factor=80, image.alpha=0.3)
dev.off()


## differential analysis

Idents(OBJ.merge)=OBJ.merge@meta.data$ImmuneLabel

## immune vs surround Immune
# group 1 over group 2, group 2 is the reference group
DEres <- FindMarkers(OBJ.merge,
                     ident.1 = "Immune",
                     ident.2 = "SurroundImmune",
                     logfc.threshold=0.05)

dim(DEres) # 16055     5

out=data.frame(Gene=rownames(DEres), DEres)
write.table(out, file="DE.immune.txt", row.names = F, col.names = T, sep="\t", quote=F)



##  surround Immune vs non immune 
# group 1 over group 2, group 2 is the reference group
DEres <- FindMarkers(OBJ.merge,
                     ident.1 = "SurroundImmune",
                     ident.2 = "NonImmune",
                     logfc.threshold=0.05)

dim(DEres) # 15732     5

out=data.frame(Gene=rownames(DEres), DEres)
write.table(out, file="DE.Surround.vs.nonimmune.txt", row.names = F, col.names = T, sep="\t", quote=F)






## only on J259X
Idents(OBJ.merge)=paste0(OBJ.merge@meta.data$library,"__", OBJ.merge@meta.data$ImmuneLabel)

# group 1 over group 2, group 2 is the reference group
DEres <- FindMarkers(OBJ.merge,
                     ident.1 = "J259X__Immune",
                     ident.2 = "J259X__SurroundImmune",
                     logfc.threshold=0.05)

dim(DEres) # 16280     5

out=data.frame(Gene=rownames(DEres), DEres)
write.table(out, file="DE.immune.J259X.txt", row.names = F, col.names = T, sep="\t", quote=F)



############# define tumor vs non-tumor area

#SpatialDimPlot(OBJ.merge, interactive = TRUE, pt.size.factor =60)

# all J259X are Tumor
TumorLabel=ifelse(OBJ.merge@meta.data$library%in%c("J259X","J259B"),"TUMOR","ADJ")
names(TumorLabel)=rownames(OBJ.merge@meta.data)

# J259B, tissue=1
posFile=OBJ.merge@images$J259B@coordinates

## select non-tumor cells
selectCell=rownames(posFile)[c(which(posFile$row==0 & posFile$col>=88),
                               which(posFile$row==1 & posFile$col>=91),
                               which(posFile$row==2 & posFile$col>=92),
                               which(posFile$row==3 & posFile$col>=93),
                               which(posFile$row==4 & posFile$col>=94),
                               which(posFile$row==5 & posFile$col>=95),
                               which(posFile$row==6 & posFile$col>=98),
                               which(posFile$row==7 & posFile$col>=105),
                               which(posFile$row==8 & posFile$col>=106),
                               which(posFile$row==9 & posFile$col>=107),
                               which(posFile$row==10 & posFile$col>=108),
                               which(posFile$row==11 & posFile$col>=109),
                               which(posFile$row==12 & posFile$col>=110),
                               which(posFile$row==13 & posFile$col>=111),
                               which(posFile$row==14 & posFile$col>=112),
                               which(posFile$row==15 & posFile$col>=113),
                               which(posFile$row==16 & posFile$col>=114),
                               which(posFile$row==17 & posFile$col>=115), 
                               which(posFile$row==18 & posFile$col>=116),
                               which(posFile$row==19 & posFile$col>=117))]

SpatialDimPlot(OBJ.merge,  pt.size.factor =60, images="J259B", cells.highlight=selectCell)


TumorLabel[selectCell]="ADJ"
table(TumorLabel)
# ADJ TUMOR 
# 238  6082 

OBJ.merge@meta.data=cbind(OBJ.merge@meta.data, TumorLabel=TumorLabel)
table(OBJ.merge@meta.data$TumorLabel)
# ADJ TUMOR 
# 238  6082 

table(OBJ.merge@meta.data$TumorLabel, OBJ.merge@meta.data$ImmuneLabel)
# Immune NonImmune SurroundImmune
# ADJ      127        27             84
# TUMOR    162      5334            586


table(paste0(OBJ.merge@meta.data$library,"__", OBJ.merge@meta.data$TumorLabel), OBJ.merge@meta.data$ImmuneLabel)
#               Immune NonImmune SurroundImmune
# J259B__ADJ      127        27             84
# J259B__TUMOR     38      1992            129
# J259X__TUMOR    124      3342            457


pdf("SpatialDimPlot.tumor.pdf", width=12, height=5)
SpatialDimPlot(OBJ.merge, group.by = "TumorLabel", pt.size.factor =60)
dev.off()

saveRDS(file = "OBJ.merge.immune.RDS", OBJ.merge)

# OBJ.merge=readRDS("OBJ.merge.immune.RDS")


pdf(paste0("DimPlot.ImmuneLabel.pdf"), width=6, height=5)
DimPlot(OBJ.merge, reduction = "umap", group.by = c("ImmuneLabel")) +
  scale_color_manual(values=c("Immune"="red", "SurroundImmune"="blue", "NonImmune"="lightgrey"))
dev.off()

pdf(paste0("DimPlot.TumorLabel.pdf"), width=6, height=5)
DimPlot(OBJ.merge, reduction = "umap", group.by = c("TumorLabel"))
dev.off()



######### DE tumor surrounding vs non-immune
Idents(OBJ.merge)=paste0(OBJ.merge@meta.data$TumorLabel,"__", OBJ.merge@meta.data$ImmuneLabel)
table(Idents(OBJ.merge))

# group 1 over group 2, group 2 is the reference group
DEres <- FindMarkers(OBJ.merge,
                     ident.1 = "TUMOR__SurroundImmune",
                     ident.2 = "TUMOR__NonImmune",
                     logfc.threshold=0.05)

dim(DEres) # 15692     5

length(which(DEres$avg_log2FC>0 & DEres$p_val<0.01)) # 773
length(which(DEres$avg_log2FC<0 & DEres$p_val<0.01)) # 429

out=data.frame(Gene=rownames(DEres), DEres)
write.table(out, file="DE.tumor_surround.vs.tumor_nonImmune.txt", row.names = F, col.names = T, sep="\t", quote=F)


######### DE adj surrounding vs non-immune
Idents(OBJ.merge)=paste0(OBJ.merge@meta.data$TumorLabel,"__", OBJ.merge@meta.data$ImmuneLabel)
table(Idents(OBJ.merge))

# group 1 over group 2, group 2 is the reference group
DEres <- FindMarkers(OBJ.merge,
                     ident.1 = "ADJ__SurroundImmune",
                     ident.2 = "ADJ__NonImmune",
                     logfc.threshold=0.05)

dim(DEres) # 16270     5

length(which(DEres$avg_log2FC>0 & DEres$p_val<0.01)) # 67
length(which(DEres$avg_log2FC<0 & DEres$p_val<0.01)) # 357

out=data.frame(Gene=rownames(DEres), DEres)
write.table(out, file="DE.adj_surround.vs.adj_nonImmune.txt", row.names = F, col.names = T, sep="\t", quote=F)



######### DE Adj surrounding immune vs tumor surrounding immune 
Idents(OBJ.merge)=paste0(OBJ.merge@meta.data$TumorLabel,"__", OBJ.merge@meta.data$ImmuneLabel)
table(Idents(OBJ.merge))

# group 1 over group 2, group 2 is the reference group
DEres <- FindMarkers(OBJ.merge,
                     ident.1 = "TUMOR__SurroundImmune",
                     ident.2 = "ADJ__SurroundImmune",
                     logfc.threshold=0.05)

dim(DEres) # 16843     5

length(which(DEres$avg_log2FC>0 & DEres$p_val<0.000001)) # 103
length(which(DEres$avg_log2FC<0 & DEres$p_val<0.000001)) # 7185

out=data.frame(Gene=rownames(DEres), DEres)
write.table(out, file="DE.tumor_surround.vs.adj_surround.txt", row.names = F, col.names = T, sep="\t", quote=F)




#################### GSVA analysis, on selected pathways
library(GSVA)
library(limma)
library(openxlsx)

rm(list=ls())

options(future.globals.maxSize= 891289600000)

SAMPLENAME=c("J259B","J259X")

OBJ.merge=readRDS("OBJ.merge.immune.RDS")

#### prepare data
geneExp=as.matrix(GetAssayData(object = OBJ.merge, slot = "data"))
dim(geneExp) # gene by cell 18085  6320
rownames(geneExp)=toupper(rownames(geneExp))
geneExp[1:5,1:5]

####### read in the pathway file and format
pathFile=read.xlsx("selected_pathway.xlsx")

PathLib=strsplit(pathFile$Genes, split=",")
names(PathLib)=pathFile$Pathway

for(i in 1:nrow(pathFile)){
  PathLib[[i]]=intersect(PathLib[[i]], rownames(geneExp))
}

### GSVA

design <- cbind(sampleGroup1=1, sampleGroup2vs1=ifelse(OBJ.merge@meta.data$TumorLabel=="TUMOR", 1, 0))

## fit linear model
fit <- lmFit(geneExp, design)

## estimate moderated t-statistics
fit <- eBayes(fit)

## genes in set1 are differentially expressed
topTable(fit, coef="sampleGroup2vs1")

## build GSVA parameter object
gsvapar <- gsvaParam(geneExp, PathLib, maxDiff=TRUE)

## estimate GSVA enrichment scores for the three sets
gsva_es <- gsva(gsvapar)
dim(gsva_es)

## output table
out=data.frame(Pathway=rownames(gsva_es),
               Gene=sapply(PathLib[rownames(gsva_es)], paste, collapse=", "),
               gsva_es)
#colnames(out)[2:10]=paste0("EnrichmentScore.",colnames(out)[2:10])
write.xlsx(out, file="GSVA.selectPath.xlsx", overwrite = T)

# out=read.xlsx("GSVA.selectPath.xlsx")

range(as.matrix(out[,3:ncol(out)])) # -0.8100450  0.6418646

## spatial plot

for(pathInd in 1:nrow(out)){
  print(out[pathInd, "Pathway"])
  
  pathwayName=out[pathInd, "Pathway"]
  pathwayName=gsub(" / ", " or ",pathwayName)
  pathwayName=gsub("/"," or ", pathwayName)
  
  OBJ.merge@meta.data$Enrichment_Score=as.numeric(out[pathInd,3:ncol(out)])
  
  pdf(paste0("AVITI.PathwayEnrichmentScore.", pathwayName, ".pdf"), width=10, height=5)
  print(SpatialFeaturePlot( OBJ.merge, features = "Enrichment_Score", 
                      pt.size.factor = 80, image.alpha = 0.3 ) )
  dev.off()
  
}


