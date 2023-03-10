library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(org.Mm.eg.db)
library(data.table)
library(harmony)
source('/data1/02.private/dengyj/analysis/mycode/DEG.Seurat.R')
gogenes <- unique(select(org.Mm.eg.db, keys = c("GO:0007049"), 
                         columns = c("SYMBOL"), keytype = "GOALL")$SYMBOL)
load('/data1/02.private/dengyj/analysis/database/mouse_cc_genes.rds')
load('/data1/02.private/dengyj/analysis/Hif2a/doublets/hif_step_1.rds')
library(openxlsx)
library(SCopeLoomR)
library(clusterProfiler)
library(ComplexHeatmap)

Idents_proj <- function(
    obj1, 
    obj2
){
    new_identity <- rbind(data.frame(b = as.character(Idents(obj1)), row.names = names(Idents(obj1))),
                    data.frame(b = as.character(Idents(obj2)), 
                    row.names = names(Idents(obj2)))[setdiff(colnames(obj2), colnames(obj1)),,drop=F])
    obj2 <- AddMetaData(obj2, metadata = new_identity, col.name ='new_identity')            
    Idents(obj2) <- obj2$new_identity
    return(obj2)
}

hif_step_1 <- readRDS('hif_step_1.rds')
hif_step_1$identity <- as.character(Idents(hif_step_1))
hif_step_1$source <- ''
hif_step_1$source[grep('^K', hif_step_1$batch)] <- 'KO'
hif_step_1$source[grep('^W', hif_step_1$batch)] <- 'WT'
hif_step_1$source[grep('^I', hif_step_1$batch)] <- 'Inh'
hif_step_1$source[grep('^O', hif_step_1$batch)] <- 'Oxy'

hif_step_1[["percent.mt"]] <- PercentageFeatureSet(hif_step_1, pattern = "mt-")
hif_step_1 <- hif_step_1[, hif_step_1$percent.mt < 10]
hif_step_1 <- hif_step_1[, !hif_step_1$batch %in% c('K2-10X5', 'K5-10X5')]
hif_step_1 <- AddMetaData(hif_step_1, metadata = orig_idents, col.name = 'orig_idents')
hif_step_1 <- NormalizeData(hif_step_1, verbose = T)
hif_step_1 <- FindVariableFeatures(hif_step_1, selection.method = "vst", nfeatures = 2000)
hif_step_1 <- CellCycleScoring(hif_step_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
hif_step_1$CC.Difference <- hif_step_1$S.Score - hif_step_1$G2M.Score
hif_step_1 <- ScaleData(hif_step_1,vars.to.regress = c('S.Score','G2M.Score'))
hif_step_1 <- RunPCA(hif_step_1, npcs = 50)
hif_step_1 <- RunHarmony(hif_step_1, 'batch')
options(repr.plot.width=6, repr.plot.height=5)
ElbowPlot(hif_step_1, ndims = 50, reduction = 'harmony')
hif_step_1 <- RunUMAP(hif_step_1, reduction = "harmony", dims = 1:30)
hif_step_1 <- FindNeighbors(hif_step_1, reduction = "harmony", dims = 1:30)
hif_step_1 <- RunTSNE(hif_step_1, reduction = "harmony", dims = 1:30)
hif_step_1 <- FindClusters(hif_step_1, resolution = 0.91)

options(repr.plot.width=6, repr.plot.height=5)
#DimPlot(hif_step_1, reduction = 'tsne',label = F) + scale_colour_manual(values = mycolors)
DimPlot(hif_step_1, reduction = 'tsne',label = T) + scale_colour_manual(values = mycolors)
DimPlot(hif_step_1, reduction = 'tsne',group.by = 'batch')+
scale_color_manual(values = mycolors)
#DimPlot(hif_step_1, reduction = 'tsne',group.by = 'orig_idents')
DimPlot(hif_step_1, reduction = 'tsne',group.by = 'Phase')
options(repr.plot.width=12, repr.plot.height=12)

Idents(hif_step_1) <- hif_step_1$RNA_snn_res.0.91
hif_step_1 <- RenameIdents(hif_step_1, '0' = 'Th2','1' ='T_memory', '2' ='Th17',  '3' = 'Treg', 
                           '4' = 'Alveolar_Mac',
                          '5' = 'IFN_T', '6' = 'Th1', '7' = 'Cd8_T', 
                          '9' = 'Trans_mac', '10' = 'T_naive', '11' = 'Neutrophil', '12' = 'Bcell', 
                           '13' ='Inflammatory_Mac','14' = 'Mki67_T',
                           '15' = 'Eosinophil','16' = 'Interstitial_mac',
                           '18' = 'Epi', 
                          '19' ='cDC', '20' = 'γδT', '21' ='Mast', '22' = "pDC",  
                           '23' = 'Epi',
                           '24' = 'Bcell')

hif_step_2 <- hif_step_1[, !Idents(hif_step_1) %in% c('8', '17')] ##去除质量不好的群
hif_step_2 <- NormalizeData(hif_step_2, verbose = T)
hif_step_2 <- FindVariableFeatures(hif_step_2, selection.method = "vst", nfeatures = 2000)
hif_step_2 <- ScaleData(hif_step_2,vars.to.regress = c('S.Score','G2M.Score'))
hif_step_2 <- RunPCA(hif_step_2, npcs = 50)
hif_step_2 <- RunHarmony(hif_step_2, 'batch')
options(repr.plot.width=6, repr.plot.height=5)
ElbowPlot(hif_step_2, ndims = 50, reduction = 'harmony')
hif_step_2 <- RunTSNE(hif_step_2, reduction = "harmony", dims = 1:30)
hif_step_2 <- RunUMAP(hif_step_2, reduction = "harmony", dims = 1:30)
hif_step_2$identity <- as.character(Idents(hif_step_2))
hif_step_2$identity_lv2 <- hif_step_2$identity
hif_step_2$identity_lv2[hif_step_2$identity_lv2 %in% 
                        c('T_naive','T_memory','Th17','Th2','Th1',
                          'Treg','IFN_T','γδT','Cd8_T','Mki67_T')] <- 'Tcell'
Idents(hif_step_2) <- hif_step_2$identity_lv2

levels(hif_step_2) <- c('Tcell', 'Bcell', 'Neutrophil', 'Eosinophil', 'Interstitial_mac', 'Trans_mac', 
                        'Alveolar_Mac', 'Inflammatory_Mac', 'cDC', 'pDC', 'Mast', 'Epi')

########进行Hif2a 的分析
Hif2a_Th2 <- hif_step_2[, hif_step_2$identity =='Th2' & hif_step_2$source %in% c('WT', 'KO')]
Hif2a_Th2 <- NormalizeData(Hif2a_Th2, verbose = T)
Hif2a_Th2 <- FindVariableFeatures(Hif2a_Th2, selection.method = "vst", nfeatures = 2000)
Hif2a_Th2 <- ScaleData(Hif2a_Th2,vars.to.regress = c('nCount_RNA'))
Hif2a_Th2 <- RunPCA(Hif2a_Th2, npcs = 50)
Hif2a_Th2 <- RunHarmony(Hif2a_Th2, 'batch')
Hif2a_Th2 <- RunTSNE(Hif2a_Th2, reduction = "harmony", dims = 1:30)
Hif2a_Th2 <- RunUMAP(Hif2a_Th2, reduction = "harmony", dims = 1:30)
Hif2a_Th2 <- FindNeighbors(Hif2a_Th2, reduction = "harmony", dims = 1:30)

DefaultAssay(Hif2a_Th2) <- 'RNA'
Hif2a_Th2 <- ScaleData(Hif2a_Th2, features = rownames(Hif2a_Th2))
markers_Hif2a_Th2 <- FindAllMarkers(Hif2a_Th2, pseudocount.use = 0, logfc.threshold = log(1.5))


Hif2a_Th2_2 <- Hif2a_Th2[, !Hif2a_Th2$seurat_clusters %in% c('5')] ###去除质量不好
Hif2a_Th2_2 <- NormalizeData(Hif2a_Th2_2, verbose = T)
Hif2a_Th2_2 <- FindVariableFeatures(Hif2a_Th2_2, selection.method = "vst", nfeatures = 2000)
Hif2a_Th2_2 <- ScaleData(Hif2a_Th2_2,vars.to.regress = c('S.Score', 'G2M.Score'))
Hif2a_Th2_2 <- RunPCA(Hif2a_Th2_2, npcs = 50)
Hif2a_Th2_2 <- RunHarmony(Hif2a_Th2_2, 'batch')
options(repr.plot.width=6, repr.plot.height=5)
ElbowPlot(Hif2a_Th2_2, ndims = 50, reduction = 'harmony')
Hif2a_Th2_2 <- RunUMAP(Hif2a_Th2_2, reduction = "harmony", dims = 1:30)
Hif2a_Th2_2 <- FindNeighbors(Hif2a_Th2_2, reduction = "harmony", dims = 1:30)
Hif2a_Th2_2 <- FindClusters(Hif2a_Th2_2, resolution = 0.45)
DefaultAssay(Hif2a_Th2_2) <- 'RNA'
Hif2a_Th2_2 <- ScaleData(Hif2a_Th2_2, features = rownames(Hif2a_Th2_2))
markers_Hif2a_Th2_2 <- FindAllMarkers(Hif2a_Th2_2, pseudocount.use = 0, logfc.threshold = log(1.5))

Hif2a_Th2<-Hif2a_Th2_2 
Hif2a_Th2$class<-Hif2a_Th2$RNA_snn_res.0.45
Hif2a_Th2$finally_class<-Hif2a_Th2$class
Hif2a_Th2$finally_class[Hif2a_Th2$class=="3"]<-"0"
Hif2a_Th2$finally_class[Hif2a_Th2$class=="0"]<-"1"
Hif2a_Th2$finally_class[Hif2a_Th2$class=="1"]<-"3"
Idents(Hif2a_Th2)<-"finally_class"

#####抑制剂
Th2_inh <- HIF[, HIF$identity =='Th2' & HIF$source %in% c("Inh")]
Th2_inh <- NormalizeData(Th2_inh, verbose = T)
Th2_inh <- FindVariableFeatures(Th2_inh, selection.method = "vst", nfeatures = 2000)
Th2_inh <- ScaleData(Th2_inh,vars.to.regress = c('nCount_RNA'))
Th2_inh <- RunPCA(Th2_inh, npcs = 50)
Th2_inh <- RunHarmony(Th2_inh, 'batch')
options(repr.plot.width=6, repr.plot.height=5)
ElbowPlot(Th2_inh, ndims = 50, reduction = 'harmony')
Th2_inh <- RunTSNE(Th2_inh, reduction = "harmony", dims = 1:30)
Th2_inh <- RunUMAP(Th2_inh, reduction = "harmony", dims = 1:30)
Th2_inh <- FindNeighbors(Th2_inh, reduction = "harmony", dims = 1:30)

Th2_inh <- FindClusters(Th2_inh, resolution = 0.5)
options(repr.plot.width=6, repr.plot.height=6)
DimPlot(Th2_inh, reduction = 'umap',label = T) + 
scale_colour_manual(values = mycolors)
DimPlot(Th2_inh, reduction = 'umap',group.by = 'batch')
DimPlot(Th2_inh, reduction = 'umap',group.by = 'Phase')
DimPlot(Th2_inh, reduction = 'umap',group.by = 'source')
Th2_inh_2<-Th2_inh[,WhichCells(Th2_inh,idents = c("0","1","2","3","4","6","7"))]

inh_I3<-Th2_inh_2[,Th2_inh_2$batch%in%"I3-10X5"]
inh_I4<-Th2_inh_2[,Th2_inh_2$batch%in%"I4-10X5"]
inh_I5<-Th2_inh_2[,Th2_inh_2$batch%in%"I5-10X5"]
inh_I6<-Th2_inh_2[,Th2_inh_2$batch%in%"I6-10X5"]
inh_I7<-Th2_inh_2[,Th2_inh_2$batch%in%"I7-10X5"]
merge_I3<-merge(m,inh_I3)
merge_I4<-merge(m,inh_I4)
merge_I5<-merge(m,inh_I5)
merge_I6<-merge(m,inh_I6)
merge_I7<-merge(m,inh_I7)
merge_all<-merge(m,list(inh_I3,inh_I4,inh_I5,inh_I6,inh_I7))
pancreas.list <- SplitObject(merge_I3, split.by = "batch")
pancreas.list <- pancreas.list[c("K3-10X5", "K4-10X5", "W1-10X5","W2-10X5", "W3-10X5","W4-10X5","I3-10X5")]
for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
    pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000,
        verbose = FALSE)
}
reference.list <- pancreas.list[c("K3-10X5", "K4-10X5", "W1-10X5","W2-10X5", "W3-10X5","W4-10X5")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)


DefaultAssay(pancreas.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)

pancreas.query <- pancreas.list[["I3-10X5"]]
pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query,
    dims = 1:30)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$class,
    dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)

I3_predclass<-pancreas.query

pancreas.list <- SplitObject(merge_I4, split.by = "batch")
pancreas.list <- pancreas.list[c("K3-10X5", "K4-10X5", "W1-10X5","W2-10X5", "W3-10X5","W4-10X5","I4-10X5")]
for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
    pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000,
        verbose = FALSE)
}
reference.list <- pancreas.list[c("K3-10X5", "K4-10X5", "W1-10X5","W2-10X5", "W3-10X5","W4-10X5")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
pancreas.query <- pancreas.list[["I4-10X5"]]
pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query,
    dims = 1:30)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$class,
    dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)
I4_predclass<-pancreas.query

pancreas.list <- SplitObject(merge_I5, split.by = "batch")
pancreas.list <- pancreas.list[c("K3-10X5", "K4-10X5", "W1-10X5","W2-10X5", "W3-10X5","W4-10X5","I5-10X5")]
for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
    pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000,
        verbose = FALSE)
}
reference.list <- pancreas.list[c("K3-10X5", "K4-10X5", "W1-10X5","W2-10X5", "W3-10X5","W4-10X5")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)

pancreas.list <- SplitObject(merge_I5, split.by = "batch")
pancreas.list <- pancreas.list[c("K3-10X5", "K4-10X5", "W1-10X5","W2-10X5", "W3-10X5","W4-10X5","I5-10X5")]
for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
    pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000,
        verbose = FALSE)
}
reference.list <- pancreas.list[c("K3-10X5", "K4-10X5", "W1-10X5","W2-10X5", "W3-10X5","W4-10X5")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
pancreas.query <- pancreas.list[["I5-10X5"]]
pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query,
    dims = 1:30,k.filter = 140)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$class,
    dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)
I5_predclass<-pancreas.query

pancreas.list <- SplitObject(merge_I6, split.by = "batch")
pancreas.list <- pancreas.list[c("K3-10X5", "K4-10X5", "W1-10X5","W2-10X5", "W3-10X5","W4-10X5","I6-10X5")]
for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
    pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000,
        verbose = FALSE)
}
reference.list <- pancreas.list[c("K3-10X5", "K4-10X5", "W1-10X5","W2-10X5", "W3-10X5","W4-10X5")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
pancreas.query <- pancreas.list[["I6-10X5"]]
pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query,
    dims = 1:30)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$class,
    dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)
I6_predclass<-pancreas.query

pancreas.list <- SplitObject(merge_I7, split.by = "batch")
pancreas.list <- pancreas.list[c("K3-10X5", "K4-10X5", "W1-10X5","W2-10X5", "W3-10X5","W4-10X5","I7-10X5")]
for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
    pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000,
        verbose = FALSE)
}
reference.list <- pancreas.list[c("K3-10X5", "K4-10X5", "W1-10X5","W2-10X5", "W3-10X5","W4-10X5")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
pancreas.query <- pancreas.list[["I7-10X5"]]
pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query,
    dims = 1:30)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$class,
    dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)
I7_predclass<-pancreas.query
Th2$predicted.id<-Th2$class
Idents(Hif2a_Th2_3)<-"predicted.id"
object<-merge(Th2,list(I3_predclass,I4_predclass,I5_predclass,I6_predclass,I7_predclass))

Hif2a_Th2_3<-object
Hif2a_Th2_3 <- NormalizeData(Hif2a_Th2_3, verbose = T)
Hif2a_Th2_3 <- FindVariableFeatures(Hif2a_Th2_3, selection.method = "vst", nfeatures = 2000)
Hif2a_Th2_3 <- ScaleData(Hif2a_Th2_3,vars.to.regress = c('nCount_RNA'))
Hif2a_Th2_3 <- RunPCA(Hif2a_Th2_3, npcs = 50)
Hif2a_Th2_3 <- RunHarmony(Hif2a_Th2_3, 'batch')
Hif2a_Th2_3 <- RunTSNE(Hif2a_Th2_3, reduction = "harmony", dims = 1:30)
Hif2a_Th2_3 <- RunUMAP(Hif2a_Th2_3, reduction = "harmony", dims = 1:30)
Hif2a_Th2_3 <- FindNeighbors(Hif2a_Th2_3, reduction = "harmony", dims = 1:30)
options(repr.plot.width=6, repr.plot.height=6)
DimPlot(Hif2a_Th2_3, reduction = 'umap',label = T,label.size = 8) + 
scale_colour_manual(values = mycolors)
