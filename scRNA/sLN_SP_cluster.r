####python#####
####去除双细胞
import numpy as np
import doubletdetection
import tarfile
import os,time
from multiprocessing import Process
file=("./filtered_feature_bc_matrix/matrix.mtx.gz")
raw_counts = doubletdetection.load_mtx(file)
zero_genes = (np.sum(raw_counts, axis=0) == 0).A.ravel()
raw_counts = raw_counts[:, ~zero_genes]
x = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=False, standard_scaling=True)
doublets = x.fit(raw_counts).predict(p_thresh=1e-16, voter_thresh=0.5)
np.savetxt("filtered_feature_bc_matrix.txt", doublets, fmt="%d", delimiter=",")

#####分群定义群
library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(org.Mm.eg.db)
library(data.table)
library(harmony)
#source('/data1/02.private/dengyj/analysis/mycode/DEG.Seurat.R')
load("/data1/02.private/guanpb/data/cellranger_DB/mouse_cc_genes_for_cellranger_2.rds")
#load('/data1/02.private/dengyj/analysis/Hif2a/doublets/hif_step_1.rds')
library(openxlsx)
library(SCopeLoomR)
library(clusterProfiler)
library(ComplexHeatmap)



dir<-paste0("./filtered_feature_bc_matrix")
data<-Read10X(data.dir =dir)
data_rna<-data[["Gene Expression"]]
data_adt<-data[["Antibody Capture"]]
xlist<-fread(paste0("./KO_SP_filtered_feature_bc_matrix.txt"))
xnames<-colnames(data_rna)[which(xlist$V1 > 0)]
data_rna<-data_rna[,setdiff(colnames(data_rna),xnames)]
x<-CreateSeuratObject(counts=data_rna,project ="SP_KO" ,min.cells = 3,min.features = 200)
data_adt<-data_adt[,colnames(x)]
x[["ADT"]]<-CreateAssayObject(counts = data_adt)
x[["percent.mt"]]<-PercentageFeatureSet(x,pattern = "mt-")
saveRDS(x,file="./SP_KO.rds")



dir<-paste0("./filtered_feature_bc_matrix")
data<-Read10X(data.dir =dir)
data_rna<-data[["Gene Expression"]]
data_adt<-data[["Antibody Capture"]]
xlist<-fread(paste0("./WT_SP_filtered_feature_bc_matrix.txt"))
xnames<-colnames(data_rna)[which(xlist$V1 > 0)]
data_rna<-data_rna[,setdiff(colnames(data_rna),xnames)]
x<-CreateSeuratObject(counts=data_rna,project ="SP_WT" ,min.cells = 3,min.features = 200)
data_adt<-data_adt[,colnames(x)]
x[["ADT"]]<-CreateAssayObject(counts = data_adt)
x[["percent.mt"]]<-PercentageFeatureSet(x,pattern = "mt-")
saveRDS(x,file="./SP_WT.rds")

SP_KO<-RenameCells(object = SP_KO,new.names = paste0(colnames(SP_KO),"_K0"))
SP_WT<-RenameCells(object = SP_WT,new.names = paste0(colnames(SP_WT),"_WT"))

load("/data1/02.private/guanpb/data/cellranger_DB/mouse_cc_genes_for_cellranger_2.rds")
SP_genotype <- NormalizeData(SP_genotype, verbose = T)
SP_genotype <- FindVariableFeatures(SP_genotype, selection.method = "vst", nfeatures = 2000)
SP_genotype <- CellCycleScoring(SP_genotype, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
SP_genotype$CC.Difference <- SP_genotype$S.Score - SP_genotype$G2M.Score
SP_genotype <- ScaleData(SP_genotype,vars.to.regress = c('S.Score','G2M.Score'))
SP_genotype <- RunPCA(SP_genotype, npcs = 50)
SP_genotype <- RunHarmony(SP_genotype, 'orig.ident')
options(repr.plot.width=6, repr.plot.height=5)
ElbowPlot(SP_genotype, ndims = 50, reduction = 'harmony')
#hif_step_1 <- RunSPNE(hif_step_1, reduction = "harmony", dims = 1:30)
SP_genotype <- RunUMAP(SP_genotype, reduction = "harmony", dims = 1:30)
SP_genotype <- FindNeighbors(SP_genotype, reduction = "harmony", dims = 1:30)
SP_genotype <- RunTSNE(SP_genotype, reduction = "harmony", dims = 1:30)
SP_genotype <- FindClusters(SP_genotype, resolution = 0.91)
#############去除高线粒体细胞 再重新分群
VlnPlot(SP_genotype,features = c("percent.mt"))
SP_genotype_2<-SP_genotype[,SP_genotype$percent.mt<=10]
SP_genotype_2 <- NormalizeData(SP_genotype_2, verbose = T)
SP_genotype_2 <- FindVariableFeatures(SP_genotype_2, selection.method = "vst", nfeatures = 2000)
SP_genotype_2 <- CellCycleScoring(SP_genotype_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
SP_genotype_2$CC.Difference <- SP_genotype_2$S.Score - SP_genotype_2$G2M.Score
SP_genotype_2 <- ScaleData(SP_genotype_2,vars.to.regress = c('S.Score','G2M.Score'))
SP_genotype_2 <- RunPCA(SP_genotype_2, npcs = 50)
SP_genotype_2 <- RunHarmony(SP_genotype_2, 'orig.ident')
options(repr.plot.width=6, repr.plot.height=5)
ElbowPlot(SP_genotype_2, ndims = 50, reduction = 'harmony')
#hif_step_1 <- RunSPNE(hif_step_1, reduction = "harmony", dims = 1:30)
SP_genotype_2 <- RunUMAP(SP_genotype_2, reduction = "harmony", dims = 1:30)
SP_genotype_2 <- FindNeighbors(SP_genotype_2, reduction = "harmony", dims = 1:30)
SP_genotype_2 <- RunTSNE(SP_genotype_2, reduction = "harmony", dims = 1:30)
SP_genotype_2 <- FindClusters(SP_genotype_2, resolution = 0.5)
saveRDS(SP_genotype_2,"SP_genotype_2.rds")



####SLN同上分析方法
dir<-paste0("。/filtered_feature_bc_matrix")
data<-Read10X(data.dir =dir)
data_rna<-data[["Gene Expression"]]
data_adt<-data[["Antibody Capture"]]
xlist<-fread(paste0("。/KO_SLN_filtered_feature_bc_matrix.txt"))
xnames<-colnames(data_rna)[which(xlist$V1 > 0)]
data_rna<-data_rna[,setdiff(colnames(data_rna),xnames)]
x<-CreateSeuratObject(counts=data_rna,project ="SLN_KO" ,min.cells = 3,min.features = 200)
data_adt<-data_adt[,colnames(x)]
x[["ADT"]]<-CreateAssayObject(counts = data_adt)
x[["percent.mt"]]<-PercentageFeatureSet(x,pattern = "mt-")
saveRDS(x,file="./SLN_KO.rds")

dir<-paste0("./filtered_feature_bc_matrix")
data<-Read10X(data.dir =dir)
data_rna<-data[["Gene Expression"]]
data_adt<-data[["Antibody Capture"]]
xlist<-fread(paste0("./WT_SLN_filtered_feature_bc_matrix.txt"))
xnames<-colnames(data_rna)[which(xlist$V1 > 0)]
data_rna<-data_rna[,setdiff(colnames(data_rna),xnames)]
x<-CreateSeuratObject(counts=data_rna,project ="SLN_WT" ,min.cells = 3,min.features = 200)
data_adt<-data_adt[,colnames(x)]
x[["ADT"]]<-CreateAssayObject(counts = data_adt)
x[["percent.mt"]]<-PercentageFeatureSet(x,pattern = "mt-")
saveRDS(x,file="./SLN_WT.rds")

SLN_KO<-RenameCells(object = SLN_KO,new.names = paste0(colnames(SLN_KO),"_K0"))
SLN_WT<-RenameCells(object = SLN_WT,new.names = paste0(colnames(SLN_WT),"_WT"))
load("/data1/02.private/guanpb/data/cellranger_DB/mouse_cc_genes_for_cellranger_2.rds")
SLN_genotype <- NormalizeData(SLN_genotype, verbose = T)
SLN_genotype <- FindVariableFeatures(SLN_genotype, selection.method = "vst", nfeatures = 2000)
SLN_genotype <- CellCycleScoring(SLN_genotype, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
SLN_genotype$CC.Difference <- SLN_genotype$S.Score - SLN_genotype$G2M.Score
SLN_genotype <- ScaleData(SLN_genotype,vars.to.regress = c('S.Score','G2M.Score'))
SLN_genotype <- RunPCA(SLN_genotype, npcs = 50)
SLN_genotype <- RunHarmony(SLN_genotype, 'orig.ident')
options(repr.plot.width=6, repr.plot.height=5)
ElbowPlot(SLN_genotype, ndims = 50, reduction = 'harmony')
SLN_genotype <- RunUMAP(SLN_genotype, reduction = "harmony", dims = 1:30)
SLN_genotype <- FindNeighbors(SLN_genotype, reduction = "harmony", dims = 1:30)
SLN_genotype <- RunTSNE(SLN_genotype, reduction = "harmony", dims = 1:30)
SLN_genotype <- FindClusters(SLN_genotype, resolution = 0.91)
#############去除高线粒体细胞 再重新分群
VlnPlot(SLN_genotype,features = c("percent.mt"))
SLN_genotype_2<-SLN_genotype[,SLN_genotype$percent.mt<=10]

SLN_genotype_2<- NormalizeData(SLN_genotype_2, verbose = T)
SLN_genotype_2 <- FindVariableFeatures(SLN_genotype_2, selection.method = "vst", nfeatures = 2000)
SLN_genotype_2 <- CellCycleScoring(SLN_genotype_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
SLN_genotype_2$CC.Difference <- SLN_genotype_2$S.Score - SLN_genotype_2$G2M.Score
SLN_genotype_2 <- ScaleData(SLN_genotype_2,vars.to.regress = c('S.Score','G2M.Score'))
SLN_genotype_2 <- RunPCA(SLN_genotype_2, npcs = 50)
SLN_genotype_2 <- RunHarmony(SLN_genotype_2, 'orig.ident')
options(repr.plot.width=6, repr.plot.height=5)
ElbowPlot(SLN_genotype_2, ndims = 50, reduction = 'harmony')
SLN_genotype_2 <- RunUMAP(SLN_genotype_2, reduction = "harmony", dims = 1:30)
SLN_genotype_2 <- FindNeighbors(SLN_genotype_2, reduction = "harmony", dims = 1:30)
SLN_genotype_2 <- RunTSNE(SLN_genotype_2, reduction = "harmony", dims = 1:30)
SLN_genotype_2 <- FindClusters(SLN_genotype_2, resolution = 0.5)
saveRDS(SLN_genotype_2,"SLN_genotype_2.rds")

###整合分析
SP<-readRDS("SP_genotype_2.rds")
SLN<-readRDS("SLN_genotype_2.rds")
HIF<-readRDS("HIF_ft_0917.rds")
SP<-RenameCells(object = SP,new.names = paste0(colnames(SP),"_SP"))
SLN<-RenameCells(object = SLN,new.names = paste0(colnames(SLN),"_SLN"))
HIF_BALF<-RenameCells(object = HIF_BALF,new.names = paste0(colnames(HIF_BALF),"_BF"))
SP$tissue<-"Sp"
SLN$tissue<-"SLN"
HIF_BALF$tissue<-"BALF"
object<-merge(HIF_BALF,list(SP,SLN))

load("/data1/02.private/guanpb/data/cellranger_DB/mouse_cc_genes_for_cellranger_2.rds")
object <- NormalizeData(object, verbose = T)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
object <- CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
object$CC.Difference <- object$S.Score - object$G2M.Score
object <- ScaleData(object,vars.to.regress = c('S.Score','G2M.Score'))
object <- RunPCA(object, npcs = 50)
object <- RunHarmony(object, 'orig.ident')
options(repr.plot.width=6, repr.plot.height=5)
ElbowPlot(object, ndims = 50, reduction = 'harmony')

object <- RunUMAP(object, reduction = "harmony", dims = 1:25)
object <- FindNeighbors(object, reduction = "harmony", dims = 1:25)
object <- RunTSNE(object, reduction = "harmony", dims = 1:25)
object <- FindClusters(object, resolution = 0.95)
options(repr.plot.width=16, repr.plot.height=10)
DimPlot(object,label = T,group.by = "identity",label.size = 7)
options(repr.plot.width=12, repr.plot.height=12)
DimPlot(object,label = T)

options(future.globals.maxSize = 10 * 1024^3)
plan('multisession', workers = 10)
DefaultAssay(object) <- 'RNA'
object <- ScaleData(object, features = rownames(object))
markers_object <- FindAllMarkers(object, pseudocount.use = 0, logfc.threshold = log(1.5))
saveRDS(markers_object,"markers_object.rds")
plan('sequential')

##################定义群
C18<-object[,WhichCells(object,idents = 18)]
C18
load("/data1/02.private/guanpb/data/cellranger_DB/mouse_cc_genes_for_cellranger_2.rds")
C18 <- NormalizeData(C18, verbose = T)
C18 <- FindVariableFeatures(C18, selection.method = "vst", nfeatures = 2000)
C18 <- CellCycleScoring(C18, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
C18$CC.Difference <- C18$S.Score - C18$G2M.Score
C18 <- ScaleData(C18,vars.to.regress = c('S.Score','G2M.Score'))
C18 <- RunPCA(C18, npcs = 50)
C18 <- RunHarmony(C18, 'orig.ident')
options(repr.plot.width=6, repr.plot.height=5)
ElbowPlot(C18, ndims = 50, reduction = 'harmony')
#hif_step_1 <- RunSPNE(hif_step_1, reduction = "harmony", dims = 1:30)
C18 <- RunUMAP(C18, reduction = "harmony", dims = 1:10)
C18 <- FindNeighbors(C18, reduction = "harmony", dims = 1:10)
C18 <- RunTSNE(C18, reduction = "harmony", dims = 1:10)
C18 <- FindClusters(C18, resolution = 0.1)
Mast<-WhichCells(C18,idents = c("5"))
Eosinphil<-WhichCells(C18,idents = c("1","3"))
Neutrophil<-WhichCells(C18,idents=c("2","4","6","7"))
C18_ftcell<-WhichCells(C18,idents=c("0"))

C69<-object[,WhichCells(object,idents = c("6","9"))]
load("/data1/02.private/guanpb/data/cellranger_DB/mouse_cc_genes_for_cellranger_2.rds")
C69 <- NormalizeData(C69, verbose = T)
C69 <- FindVariableFeatures(C69, selection.method = "vst", nfeatures = 2000)
C69 <- CellCycleScoring(C69, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
C69$CC.Difference <- C69$S.Score - C69$G2M.Score
C69 <- ScaleData(C69,vars.to.regress = c('S.Score','G2M.Score'))
C69 <- RunPCA(C69, npcs = 50)
C69 <- RunHarmony(C69, 'orig.ident')
options(repr.plot.width=6, repr.plot.height=5)
ElbowPlot(C69, ndims = 50, reduction = 'harmony')
#hif_step_1 <- RunSPNE(hif_step_1, reduction = "harmony", dims = 1:30)
C69 <- RunUMAP(C69, reduction = "harmony", dims = 1:20)
C69 <- FindNeighbors(C69, reduction = "harmony", dims = 1:20)
C69 <- RunTSNE(C69, reduction = "harmony", dims = 1:20)
C69 <- FindClusters(C69, resolution = 1.2)
C69_ftcell<-WhichCells(C69,idents=c("10","11"))

C16<-object[,WhichCells(object,idents = c("16"))]
load("/data1/02.private/guanpb/data/cellranger_DB/mouse_cc_genes_for_cellranger_2.rds")
C16 <- NormalizeData(C16, verbose = T)
C16 <- FindVariableFeatures(C16, selection.method = "vst", nfeatures = 2000)
C16 <- CellCycleScoring(C16, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
C16$CC.Difference <- C16$S.Score - C16$G2M.Score
C16 <- ScaleData(C16,vars.to.regress = c('S.Score','G2M.Score'))
C16 <- RunPCA(C16, npcs = 50)
C16 <- RunHarmony(C16, 'orig.ident')
options(repr.plot.width=6, repr.plot.height=5)
ElbowPlot(C16, ndims = 50, reduction = 'harmony')
C16 <- RunUMAP(C16, reduction = "harmony", dims = 1:10)
C16 <- FindNeighbors(C16, reduction = "harmony", dims = 1:10)
C16 <- RunTSNE(C16, reduction = "harmony", dims = 1:10)
C16 <- FindClusters(C16, resolution = 0.6)


Cd8_T<-WhichCells(C69,idents=c("12","9","4","16"))
Th1<-WhichCells(C69,idents = c("6","7","2","13","14"))
gdT<-WhichCells(C69,idents=c("15"))
Th17<-WhichCells(C69,idents = c("5","18","1","0","3","17","8"))
Alveolar_Mac<-WhichCells(C16,idents=c("2","3","9"))
Interstitial_Mac<-WhichCells(C16,idents=c("0","7"))
Trans_Mac<-WhichCells(C16,idents=c("4","8","6","1","5"))



###############
object$raw_idty<-object$RNA_snn_res.0.95
object$new_idty[colnames(object)%in%WhichCells(object,idents=c("24"))]<-"Epi"
object$new_idty[colnames(object)%in%WhichCells(object,idents=c("0","1","20","22"))]<-"T_naive"
object$new_idty[colnames(object)%in%WhichCells(object,idents=c("2","4","11"))]<-"T_memory"
object$new_idty[colnames(object)%in%WhichCells(object,idents=c("17","15"))]<-"Proliferation"
object$new_idty[colnames(object)%in%WhichCells(object,idents=c("10","25"))]<-"Plasma cell"
object$new_idty[colnames(object)%in%WhichCells(object,idents=c("13"))]<-"B cell"
object$new_idty[colnames(object)%in%WhichCells(object,idents=c("7","8"))]<-"Treg_1"
object$new_idty[colnames(object)%in%WhichCells(object,idents=c("3"))]<-"Th2"
object$new_idty[colnames(object)%in%WhichCells(object,idents=c("5"))]<-"Tfh"
object$new_idty[colnames(object)%in%WhichCells(object,idents=c("14"))]<-"Treg_2"
object$new_idty[colnames(object)%in%WhichCells(object,idents=c("12"))]<-"IFN_T"
object$new_idty[colnames(object)%in%WhichCells(object,idents=c("19"))]<-"Inflammatory_Mac"
object$new_idty[colnames(object)%in%Cd8_T]<-"Cd8_T"
object$new_idty[colnames(object)%in%Th1]<-"Th1"
object$new_idty[colnames(object)%in%gdT]<-"γδT"
object$new_idty[colnames(object)%in%Th17]<-"Th17"
object$new_idty[colnames(object)%in%WhichCells(object,idents=c("23"))]<-"pDC"
object$new_idty[colnames(object)%in%WhichCells(object,idents=c("21"))]<-"cDC"
object$new_idty[colnames(object)%in%Alveolar_Mac]<-"Alveolar_Mac"
object$new_idty[colnames(object)%in%Interstitial_Mac]<-"Interstitial_Mac"
object$new_idty[colnames(object)%in%Trans_Mac]<-"Trans_Mac"
object$new_idty[colnames(object)%in%c(C18_ftcell,C69_ftcell)]<-"ft_cell"
object$new_idty[colnames(object)%in%Mast]<-"Mast"
object$new_idty[colnames(object)%in%Eosinphil]<-"Eosinophil"
object$new_idty[colnames(object)%in%Neutrophil]<-"Neutrophil"




















