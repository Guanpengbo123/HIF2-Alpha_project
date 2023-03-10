####去双细胞同scRNA-seq处理方式############

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
doublets_list <- read.table("./Ts_SM_KO_filtered_feature_bc_matrix.txt")

iter_exp<- Read10X("./filtered_feature_bc_matrix")#[, doublets_list$V1 ==0]

iter_rna <- iter_exp$`Gene Expression`[, doublets_list$V1 ==0]

iter_object <- CreateSeuratObject(iter_rna, project ="Ts_ab_KO" , min.cells = 3, min.features = 200)

iter_object[['ADT']] <- CreateAssayObject(counts = iter_exp$`Antibody Capture`[, colnames(iter_object)])

iter_object  <- RenameCells(iter_object, new.names = paste(colnames(iter_object), "Ts_ab_KO",sep = '_'))

iter_object$batch <- "Ts_ab_KO" 

saveRDS(iter_object,file="./Ts_ab_KO.rds")

doublets_list <- read.table("./Ts_SM_WT_filtered_feature_bc_matrix.txt")

iter_exp<- Read10X("./filtered_feature_bc_matrix")#[, doublets_list$V1 ==0]

iter_rna <- iter_exp$`Gene Expression`[, doublets_list$V1 ==0]

iter_object <- CreateSeuratObject(iter_rna, project ="Ts_ab_WT" , min.cells = 3, min.features = 200)

iter_object[['ADT']] <- CreateAssayObject(counts = iter_exp$`Antibody Capture`[, colnames(iter_object)])

iter_object  <- RenameCells(iter_object, new.names = paste(colnames(iter_object), "Ts_ab_WT",sep = '_'))

iter_object$batch <- "Ts_ab_WT" 

Ts_ab_genotype<-merge(Ts_ab_KO,Ts_ab_WT)
Ts_ab_genotype<-merge(Ts_ab_KO,Ts_ab_WT)
Ts_ab_genotype[["percent.Rpsl"]]<-PercentageFeatureSet(Ts_ab_genotype,pattern = "Rp")
Ts_ab_genotype <- NormalizeData(Ts_ab_genotype, verbose = T)
Ts_ab_genotype <- FindVariableFeatures(Ts_ab_genotype, selection.method = "vst", nfeatures = 2000)
Ts_ab_genotype <- ScaleData(Ts_ab_genotype,vars.to.regress = c('nCount_RNA'))
Ts_ab_genotype <- RunPCA(Ts_ab_genotype, npcs = 50)
Ts_ab_genotype <- RunHarmony(Ts_ab_genotype, 'orig.ident')
options(repr.plot.width=6, repr.plot.height=5)
ElbowPlot(Ts_ab_genotype, ndims = 50, reduction = 'harmony')
Ts_ab_genotype <- RunTSNE(Ts_ab_genotype, reduction = "harmony", dims = 1:30)
Ts_ab_genotype <- RunUMAP(Ts_ab_genotype, reduction = "harmony", dims = 1:30)
Ts_ab_genotype <- FindNeighbors(Ts_ab_genotype, reduction = "harmony", dims = 1:30)
Ts_ab_genotype <- FindClusters(Ts_ab_genotype, resolution = 0.91)

HIF<-readRDS("/data1/02.private/other/hif_step_2.rds")
T_idents <- c('Th2','Treg','T_memory','Th17','IFN_T','Th1','Cd8_T','T_naive','Mki67_T','γδT')
T_object<-HIF[,HIF$identity%in%T_idents]
Hif_T_Ts_ab<-merge(T_object,Ts_ab_genotype)
Hif_T_Ts_ab <- NormalizeData(Hif_T_Ts_ab, verbose = T)
Hif_T_Ts_ab <- FindVariableFeatures(Hif_T_Ts_ab, selection.method = "vst", nfeatures = 2000)
Hif_T_Ts_ab <- CellCycleScoring(Hif_T_Ts_ab, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
Hif_T_Ts_ab$CC.Difference <- Hif_T_Ts_ab$S.Score - Hif_T_Ts_ab$G2M.Score
Hif_T_Ts_ab <- ScaleData(Hif_T_Ts_ab,vars.to.regress = c('S.Score','G2M.Score'))
Hif_T_Ts_ab <- RunPCA(Hif_T_Ts_ab, npcs = 50)
Hif_T_Ts_ab <- RunHarmony(Hif_T_Ts_ab, 'orig.ident')
options(repr.plot.width=6, repr.plot.height=5)
ElbowPlot(Hif_T_Ts_ab, ndims = 50, reduction = 'harmony')
#hif_step_1 <- RunTSNE(hif_step_1, reduction = "harmony", dims = 1:30)
Hif_T_Ts_ab <- RunUMAP(Hif_T_Ts_ab, reduction = "harmony", dims = 1:30)
Hif_T_Ts_ab <- FindNeighbors(Hif_T_Ts_ab, reduction = "harmony", dims = 1:30)
Hif_T_Ts_ab <- FindClusters(Hif_T_Ts_ab, resolution = 0.91)
Hif_T_Ts_ab <- RunTSNE(Hif_T_Ts_ab, reduction = "harmony", dims = 1:30)

Th2_object<-Hif_T_Ts_ab[,WhichCells(Hif_T_Ts_ab,idents =c("3","6","22","21","18","19","9"))]

Th2<-readRDS("/data1/02.private/other/Hif2a_Th2_2.rds")
Th2$class<-Th2$RNA_snn_res.0.45
Th2_2<-merge(m,Ts_ab_Th2_2)
Th2_2 <- NormalizeData(Th2_2, verbose = T)
Th2_2 <- FindVariableFeatures(Th2_2, selection.method = "vst", nfeatures = 2000)
Th2_2 <- ScaleData(Th2_2,vars.to.regress = c('nCount_RNA'))
Th2_2 <- RunPCA(Th2_2, npcs = 50)
Th2_2 <- RunHarmony(Th2_2, 'batch')
options(repr.plot.width=6, repr.plot.height=5)
ElbowPlot(Th2_2, ndims = 50, reduction = 'harmony')
Th2_2 <- RunTSNE(Th2_2, reduction = "harmony", dims = 1:30)
Th2_2 <- RunUMAP(Th2_2, reduction = "harmony", dims = 1:30)
Th2_2 <- FindNeighbors(Th2_2, reduction = "harmony", dims = 1:30)
Th2_2 <- FindClusters(Th2_2, resolution = 0.5)
options(repr.plot.width=6, repr.plot.height=6)
DimPlot(Th2_2, reduction = 'umap',label = T,group.by = "class") + 
scale_colour_manual(values = mycolors)
#DimPlot(Th2_2, reduction = 'umap',group.by = 'batch')
#DimPlot(Th2_2, reduction = 'umap',group.by = 'Phase')
#DimPlot(Th2_2, reduction = 'umap',group.by = 'source')

#Th2_2 <- RunTSNE(Th2_2, reduction = "harmony", dims = 1:30)
Th2_2 <- RunUMAP(Th2_2, reduction = "harmony", dims = 1:30)
Th2_2 <- FindNeighbors(Th2_2, reduction = "harmony", dims = 1:30)
Th2_2 <- FindClusters(Th2_2, resolution = 0.5)
options(repr.plot.width=6, repr.plot.height=6)
DimPlot(Th2_2, reduction = 'umap',label = T,group.by = "class") + 
scale_colour_manual(values = mycolors)
DimPlot(Th2_2, reduction = 'umap',label = T) + 
scale_colour_manual(values = mycolors)

Th2_2$use<-""
Th2_2$use[WhichCells(Th2_2,idents = c("4","9"))]<-"3"
Th2_2$use[WhichCells(Th2_2,idents = c("0","6","5","10"))]<-"0"
Th2_2$use[WhichCells(Th2_2,idents = c("1","8"))]<-"1"
Th2_2$use[WhichCells(Th2_2,idents = c("2"))]<-"2"
Th2_2$use[WhichCells(Th2_2,idents = c("3","7"))]<-"4"
Th2_2$use[WhichCells(Th2_2,idents = c("11"))]<-"5"

DefaultAssay(Th2_2)<-"ADT"

Th2_2 <- NormalizeData(Th2_2, verbose = T,normalization.method = 'CLR')
Th2_2 <- ScaleData(Th2_2)

Th2_2$fin_class<-Th2_2$use
Th2_2$fin_class[Th2_2$use=="0"]<-"1"
Th2_2$fin_class[Th2_2$use=="1"]<-"3"
Th2_2$fin_class[Th2_2$use=="3"]<-"0"
