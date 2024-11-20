library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(harmony)

library(data.table)

library(reshape2)
library(plyr)

library(Matrix)

library(randomForest) 
library(pROC)

library(openxlsx)

library(patchwork)


library(ggforce)

library(ComplexHeatmap)

dim(combined_CD4T)

library(VISION)

source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')


source('/data1/02.private/dengyj/analysis/thymus/plot/support/plot_helper.r')

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

load('/data1/02.private/dengyj/analysis/database/mouse_cc_genes_for_cellranger_3.rds')

GSE131935 <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/public_data/Mm_Th2_Immunity_2019/GSE131935_new_id.rds')

levels(GSE131935)

# GSE190795 <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/public_data/GSE190795//GSE190795.rds')

# levels(GSE190795)

GSE196470_CD4T <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/public_data/GSE196470//CD4T.rds')

levels(GSE196470_CD4T)

common_features <- intersect(rownames(GSE131935), rownames(GSE196470_CD4T))
#common_features <- intersect(common_features, rownames())
length(common_features)

GSE131935$Identity <- as.character(Idents(GSE131935))
GSE196470_CD4T$Identity <- as.character(Idents(GSE196470_CD4T))

combined_CD4T <- merge(GSE131935[, !Idents(GSE131935) %in% c('Myeloid', 'undefined')], 
                         
                             GSE196470_CD4T[, !grepl('Ascr',GSE196470_CD4T$batch)])

combined_CD4T[['common']] <- CreateAssayObject(counts = combined_CD4T$RNA@counts[common_features, ])

combined_CD4T$common@data <- combined_CD4T$RNA@data[rownames(combined_CD4T$common), ]

combined_CD4T$source <- ''
combined_CD4T$source[colnames(combined_CD4T) %in% colnames(GSE131935)] <- 'GSE131935'
#combined_CD4T$source[colnames(combined_CD4T) %in% colnames(GSE190795)] <- 'GSE190795'
combined_CD4T$source[colnames(combined_CD4T) %in% colnames(GSE196470_CD4T)] <- 'GSE196470'

DefaultAssay(combined_CD4T) <- 'common'
combined_CD4T <- FindVariableFeatures(combined_CD4T,selection.method = "vst", nfeatures = 2000)
#combined_CD4T <- CellCycleScoring(combined_CD4T, s.features = s.genes, g2m.features = g2m.genes)
combined_CD4T <- ScaleData(combined_CD4T)#, vars.to.regress = c('S.Score', 'G2M.Socre'))
combined_CD4T <- RunPCA(combined_CD4T, npcs = 50, verbose = T)

combined_CD4T <- RunHarmony(combined_CD4T, c('batch'),assay.use = 'common')



options(repr.plot.width=6, repr.plot.height=6)
ElbowPlot(combined_CD4T, reduction = 'harmony',ndims = 50)


#combined_CD4T <- RunTSNE(combined_CD4T, reduction = "pca", dims = setdiff(1:20, c()), verbose = F)

combined_CD4T <- RunUMAP(combined_CD4T, reduction = "harmony", dims = setdiff(1:15, c()), verbose = F)

options(repr.plot.width=10, repr.plot.height=8)
DimPlot(combined_CD4T, reduction = 'umap', label = T, group.by = 'batch')

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(combined_CD4T, reduction = 'tsne', label = T)

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(combined_CD4T, reduction = 'umap', label = T, group.by = 'source')

options(repr.plot.width=12, repr.plot.height=6)
DimPlot(combined_CD4T, reduction = 'umap', label = T, group.by = 'Identity', split.by = 'source')

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(combined_CD4T, reduction = 'umap', label = T, group.by = 'Phase')

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(combined_CD4T) <- 'RNA'
FeaturePlot(combined_CD4T, features =c ('percent.mt', 'nCount_RNA', 'nFeature_RNA'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(combined_CD4T) <- 'RNA'
FeaturePlot(combined_CD4T, features =c ('Tmem176a', 'Tmem176b', 'Rorc', 'Ccr6', 'Il17re', 'Ccr4', 'Il17a', 'Il17f'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(combined_CD4T) <- 'RNA'
FeaturePlot(combined_CD4T, features =c ('Tbx21', 'Ifng', 'Ccl5', 'Gzmk'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(combined_CD4T) <- 'RNA'
FeaturePlot(combined_CD4T, features =c ('Gata3', 'Il1rl1', 'Tcf7', 'Lef1', 'Epas1', 'Ccr7', 'Il4', 'Il5', 'Il13'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(combined_CD4T) <- 'RNA'
FeaturePlot(combined_CD4T, features =c ('Foxp3', 'Ikzf2', 'Ctla4', 'Tnfrsf8'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(combined_CD4T) <- 'RNA'
FeaturePlot(combined_CD4T, features =c ('Ahnak', 'Cd44', 'Anxa2', 'Anxa1'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(combined_CD4T) <- 'RNA'
FeaturePlot(combined_CD4T, features =c ('Sell'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(combined_CD4T) <- 'RNA'
FeaturePlot(combined_CD4T, features =c ('Mki67', 'Pcna', 'Birc5'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(combined_CD4T) <- 'RNA'
FeaturePlot(combined_CD4T, features =c ('Srm', 'Wdr75', 'C1qbp', 'Rpf2'), 
            cols = c('lightgrey', 'red'))



cluster_colors <- 
c("Th2"="#E5D2DD","T_memory"="#53A85F","Th17"="#F1BB72","Treg"="#F3B1A0",'ST2+Treg'='#859472',
  "Alveolar_Mac"="#D6E7A3",
"IFN_T"="#57C3F3","Th1"="#476D87","Cd8_T"="#E95C59","Trans_mac"="#E59CC4","T_naive"="#AB3282",
"Neutrophil"="#23452F","Bcell"="#BD956A","Inflammatory_Mac"="#8C549C","Mki67_T"="#585658",
"Eosinophil"="#9FA3A8","Interstitial_mac"="#E0D4CA","Epi"="#5F3D69","cDC"="#C5DEBA",
  "γδT"="#58A4C3",
"Mast"="#E4C755","pDC"="#F7F398","17"="#AA9A59","8"="#E63863","other" = "gray80", 
 'Tcell'='#4284c9', 
 "0"="#C61584","1"="#AAD0BA","2"="#5CA553","3"="#E78C5D","4"="#66AAE2","5"="#F9EE40")

DefaultAssay(combined_CD4T) <- 'common'
combined_CD4T <- FindNeighbors(combined_CD4T, reduction = "harmony", dims = 1:15)


combined_CD4T <- FindClusters(combined_CD4T, resolution = seq(0.1, 1, 0.05), verbose = F)

options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.1, 1, 0.05)){
    res <- paste0('common_snn_res.', i)
    p <- DimPlot(combined_CD4T, reduction = 'umap',label = T, label.size = 20/.pt, 
                 group.by = res)

    print(p)
    print(res)
}



Idents(combined_CD4T) <- combined_CD4T$common_snn_res.0.5
options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combined_CD4T, reduction = 'umap',label = T, label.size = 20/.pt)

saveRDS(colnames(combined_CD4T)[combined_CD4T$common_snn_res.0.5 == '5'], file = 'cells.rds')

C4 <- combined_CD4T[, Idents(combined_CD4T) %in% c('4')]

DefaultAssay(C4) <- 'common'

C4 <- FindNeighbors(C4, reduction = "harmony", dims = 1:10)


C4 <- FindClusters(C4,resolution = 0.1, verbose = F)


options(repr.plot.width=8, repr.plot.height=8)
DimPlot(C4, reduction = 'umap',label = T, label.size = 20/.pt)

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(C4) <- 'RNA'
FeaturePlot(C4, features =c ('Isg15', 'Ifng', 'Tbx21'), 
            cols = c('lightgrey', 'red'))

Idents(C4) <- C4$seurat_clusters
options(repr.plot.width=8, repr.plot.height=8)
C4 <- RenameIdents(C4, '0'='Th1', '1' = 'IFN_T')
DimPlot(C4, reduction = 'umap',label = T, label.size = 20/.pt)



options(repr.plot.width = 8, repr.plot.height = 8)
FeaturePlot(Th2, features =c('Klrb1', 'Ptgdr2', 'Hpgds'))


CD4_sub <- combined_CD4T[, Idents(combined_CD4T) %in% c('5')]
DefaultAssay(CD4_sub) <- 'common'

options(repr.plot.width=6, repr.plot.height=6)
ElbowPlot(CD4_sub, ndims = 50, reduction = 'harmony')

CD4_sub <- RunUMAP(CD4_sub, reduction = "harmony", dims = 1:10, verbose = F)

options(repr.plot.width=16, repr.plot.height=8)

DimPlot(CD4_sub, reduction = 'umap',label = T, label.size = 20/.pt, group.by = 'Identity', split.by = 'source')

options(repr.plot.width=8, repr.plot.height=8)

DimPlot(CD4_sub, reduction = 'umap',label = T, label.size = 20/.pt, group.by = 'batch')

options(repr.plot.width=8, repr.plot.height=8)

DimPlot(CD4_sub, reduction = 'umap',label = T, label.size = 20/.pt, group.by = 'source')

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(CD4_sub) <- 'RNA'
FeaturePlot(CD4_sub, features =c ('percent.mt', 'nCount_RNA', 'nFeature_RNA'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(CD4_sub) <- 'RNA'
FeaturePlot(CD4_sub, features =c ('Gata3', 'Il1rl1', 'Tcf7', 'Lef1', 'Epas1', 'Ccr7', 'Il4', 'Il5', 'Il13'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(CD4_sub) <- 'RNA'
FeaturePlot(CD4_sub[, CD4_sub$source =='GSE131935'], features =c ('Gata3', 'Il1rl1', 'Tcf7', 'Lef1', 'Epas1', 'Ccr7', 'Il4', 'Il5', 'Il13'), 
            cols = c('lightgrey', 'red'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(CD4_sub) <- 'RNA'
FeaturePlot(CD4_sub, features =c ( 'Tcf7', 'Lef1', 'Ccr7', 'Sell'), 
            cols = c('lightgrey', 'red'))





DefaultAssay(CD4_sub) <- 'common'
CD4_sub <- FindNeighbors(CD4_sub, reduction = "harmony", dims = 1:10)


CD4_sub <- FindClusters(CD4_sub, resolution = seq(0.1, 0.8, 0.05), verbose = F)

options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.1, 0.8, 0.05)){
    res <- paste0('common_snn_res.', i)
    p <- DimPlot(CD4_sub, reduction = 'umap',label = T, label.size = 20/.pt, 
                 group.by = res)

    print(p)
    print(res)
}

options(repr.plot.width=8, repr.plot.height=8)
Idents(CD4_sub)<- CD4_sub$`common_snn_res.0.2`
DimPlot(CD4_sub, reduction = 'umap',label = T, label.size = 20/.pt)

Idents(CD4_sub) <- CD4_sub$`common_snn_res.0.2`
CD4_sub <- RenameIdents(CD4_sub, '0' = 'ST2_hi_Th2','2' = 'low_quality','1' = 'TCM')

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(CD4_sub, reduction = 'umap',label = T, label.size = 20/.pt, pt.size = 2)+
theme(panel.border=element_rect(colour = "black", fill =NA),
          panel.background = element_rect(colour = NA, fill ='white'), 
      axis.title = element_text(size = 15),axis.text = element_text(size = 15),
      legend.title = element_text(size = 15),legend.text = element_text(size = 15),
     plot.title = element_text(size = 15, face='italic', hjust=0.5))+
guides(color = F)





options(repr.plot.width=8, repr.plot.height=8)
Idents(combined_CD4T) <- combined_CD4T$common_snn_res.0.5
DimPlot(combined_CD4T, reduction = 'umap',label = T, label.size = 20/.pt, pt.size = 1)

length(unique(combined_CD4T$`common_snn_res.0.5`))

levels(CD4_sub)

levels(C4)

levels(combined_CD4T)

Idents(combined_CD4T) <- combined_CD4T$`common_snn_res.0.5`
options(repr.plot.width=8, repr.plot.height=8)
combined_CD4T <- RenameIdents(combined_CD4T, '6' = 'Mki67+T','1' = 'ST2hi_Th2',
                              '3' = 'ST2lo_Th2',
                              '0' = 'Naive','7' = 'Th17','2' = 'Treg', '8' = 'Treg')
combined_CD4T <- Idents_proj(CD4_sub, combined_CD4T)
combined_CD4T <- Idents_proj(C4, combined_CD4T)

DimPlot(combined_CD4T, reduction = 'umap',label = T, label.size = 20/.pt, pt.size = 1)

filtered_combined_CD4T <- combined_CD4T[, !Idents(combined_CD4T) %in% c('low_quality')]

filtered_combined_CD4T$id_lv1 <- Idents(filtered_combined_CD4T)
filtered_combined_CD4T$id_lv2 <- filtered_combined_CD4T$id_lv1
levels(filtered_combined_CD4T$id_lv2)[grepl('Th2', levels(filtered_combined_CD4T$id_lv2))] <- 'Th2'

Idents(filtered_combined_CD4T) <- filtered_combined_CD4T$id_lv2

saveRDS(filtered_combined_CD4T, file = 'filtered_combined_CD4T.rds')
saveRDS(combined_CD4T, file = 'combined_CD4T.rds')