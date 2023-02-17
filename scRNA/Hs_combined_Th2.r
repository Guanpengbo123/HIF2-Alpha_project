####combined CD4T/ILC2 from human disease model(asthma and CRSwNP)
library(Seurat)
library(ggplot2)
library(future)
library(openxlsx)
library(harmony)
library(Matrix)
library(rlang)
library(plyr)
library(patchwork)
library(reshape2)

load('/data1/02.private/dengyj/analysis/database/specie_conversion/specie_conversion.rds')
source('/data1/02.private/dengyj/analysis/Hif2a/plot/support/plot_helper.r')



###to load data contain human CD4T/ILC2 from asthma and CRSwNP model
GSE179292_2 <- readRDS('../Hs_nasal_polyp_SciImm_2021/GSE179292_2.rds')
CD4_T <- readRDS('../Hs_asthma_NatMed_2019/CD4_T.rds')



options(repr.plot.width=8, repr.plot.height=8)
DimPlot(GSE179292_2, reduction = 'umap', label = T, group.by='Identity')

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(CD4_T, reduction = 'tsne', label = T, group.by='Cluster')


GSE179292_2$source <- 'GSE179292'

CD4_T$source <- 'NatMed_asthma'

common_features <- intersect(rownames(CD4_T), rownames(GSE179292_2))


length(common_features)

combinedT <- merge(GSE179292_2[common_features, !GSE179292_2$Identity %in% c('non_Lym')], 
                  #GSE175930_Tcell[common_features, GSE175930_Tcell$Identity %in% c('CD4T')], 
                      CD4_T[common_features, ])

combinedT_list <- SplitObject(combinedT, split.by = 'source')
combinedT_list <- lapply(combinedT_list, function(x){
    #x <- NormalizeData(x, verbose = F)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = F)
})

anchor_features <- SelectIntegrationFeatures(combinedT_list, 2000)

anchors <- FindIntegrationAnchors(object.list = combinedT_list, dims = 1:20, 
                                  anchor.features = anchor_features)
combinedT <- IntegrateData(anchorset = anchors, dims = 1:20)

options(future.globals.maxSize = 10 * 1024^3)
#plan('multisession', workers = 10)
DefaultAssay(combinedT) <- 'integrated'
combinedT <- ScaleData(combinedT)#, vars.to.regress = c('nCount_RNA'))
#plan('sequential')


combinedT <- RunPCA(combinedT, npcs = 50, verbose = T)

options(repr.plot.width=6, repr.plot.height=6)
ElbowPlot(combinedT, ndims = 50, reduction = 'pca')


print(combinedT[['pca']], 1:20)

DefaultAssay(combinedT) <- 'integrated'

combinedT <- RunUMAP(combinedT, reduction = "pca", dims = 1:10)

combinedT <- RunTSNE(combinedT, reduction = "pca", dims = 1:10)

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combinedT, reduction = 'umap', group.by = 'Phase')
options(repr.plot.width=12, repr.plot.height=8)
DimPlot(combinedT, reduction = 'umap', group.by = 'source')


options(repr.plot.width=16, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('nCount_RNA', 'nFeature_RNA'),ncol = 2,
            cols = c("lightgrey", 'red'))  

DefaultAssay(combinedT) <- 'integrated'
combinedT <- FindNeighbors(combinedT, reduction = "pca", dims = 1:10)
DefaultAssay(combinedT) <- 'RNA'

DefaultAssay(combinedT) <- 'integrated'
combinedT <- FindClusters(combinedT, resolution = 0.4)
DefaultAssay(combinedT) <- 'RNA'

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combinedT, label = T)

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combinedT, label = T, reduction = 'tsne')

DefaultAssay(combinedT) <- 'RNA'

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('IL13', 'GATA3', 'PTGDR2'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=16, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('S1PR1', 'KLF2'),ncol = 2,
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('NR4A2'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('CD3D', 'CD3E', 'CD3G'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('CD4', 'CD40LG', 'CD8A', 'CD8B'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('SLC4A10', 'RORC', 'NCR3', 'ZBTB16'),
            cols = c("lightgrey", 'red'))              


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('IL13', 'IL1RL1', 'GATA3', 'EPAS1',
                                   'CD200R1', 'IL5', 'IL4', 'PTGDR2', 'IL10'),
            cols = c("lightgrey", 'red'))                        

####Th2， ILC2
options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('TCN1','HPGDS','CD109','PDCD1', 'CD200', 'TOX2', 'FOXP3', 'XCL1', 'CD3E'),
            cols = c("lightgrey", 'red'))   


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('TBX21', 'IFNG', 'CCR6', 'RORC'),
            cols = c("lightgrey", 'red'))        


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('IL17A', 'IL17F', 'CCR6', 'RORC'),
            cols = c("lightgrey", 'red'))        


options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('IL23R'),
            cols = c("lightgrey", 'red'))        


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, reduction = 'umap',features = c('TCF7', 'LEF1', 'SELL', 'BATF'),
            cols = c("lightgrey", 'red'))        

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, reduction = 'umap',features = c('IL7R', 'GPR183', 'AHNAK', 'ANXA1'),
            cols = c("lightgrey", 'red')) 

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, reduction = 'umap',features = c('CD44', 'FAS', 'ANXA2', 'PASK'),
            cols = c("lightgrey", 'red')) 

options(repr.plot.width=16, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('CD44', 'FAS', 'ANXA2', 'PASK', 
                                         'IL7R', 'GPR183', 'AHNAK', 'ANXA1'),ncol=4,
            cols = c("lightgrey", 'red')) 

options(repr.plot.width=16, repr.plot.height=4)#
FeaturePlot(combinedT, reduction = 'umap',features = c('PDCD1', 'BCL6', 'TOX2', 'IL6ST'),ncol=4,
            cols = c("lightgrey", 'red')) 

options(repr.plot.width=16, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('CCR6', 'RORC'),ncol=2,
            cols = c("lightgrey", 'red')) 

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, reduction = 'umap',features = c('IFNG', 'TBX21', 'GZMK', 'CCL5'),#ncol=3,
            cols = c("lightgrey", 'red')) 


options(repr.plot.width=16, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('S1PR1','KLF2'),ncol =2,
            cols = c("lightgrey", 'red'))  


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('TCF7', 'LEF1', 'SELL', 'BATF'),
            cols = c("lightgrey", 'red'))        


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('CXCR5', 'BCL6', 'PDCD1', 'TOX2'),
            cols = c("lightgrey", 'red'))        


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('GZMB', 'GZMK', 'PRF1', 'CX3CR1'),
            cols = c("lightgrey", 'red'))        



options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('CD8A', 'CD8B', 'SLC4A10'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('CXCR6'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('CD27'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('S1PR1'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('SATB1'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combinedT, reduction = 'umap',features = c('PP1CB', 'CREM', 'LRRK2', 'YPEL5'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('IL17RC'),
            cols = c("lightgrey", 'red'))  

markers_combinedT <- FindAllMarkers(combinedT, pseudocount.use = 0, logfc.threshold = log(1.5))

openxlsx:::write.xlsx(markers_combinedT, file = 'markers_combinedT.xlsx')

DefaultAssay(combinedT) <- 'RNA'
combinedT <- ScaleData(combinedT, features = rownames(combinedT))
options(repr.plot.width=9, repr.plot.height=9)
p <- DoHeatmap(combinedT, features = unique(markers_combinedT[markers_combinedT$avg_logFC > log(2) & 
                                                    markers_combinedT$p_val_adj < 0.05 &
                                                    markers_combinedT$pct.1 > 0.2, 'gene']))
p

rownames(combinedT)[grep('IL17', rownames(combinedT))]

options(repr.plot.width=16, repr.plot.height=16)
for(i in unique(markers_combinedT$cluster)){
    genes <- markers_combinedT[markers_combinedT$cluster == i & 
                                 markers_combinedT$avg_logFC > log(1.5) &
                          markers_combinedT$p_val_adj < 0.05, 'gene'][1:16]
    genes <- genes[!is.na(genes)]
    if(length(genes) > 0){
        p <- FeaturePlot(combinedT, reduction = 'umap', cols = c('lightgrey', 'red'),
                features = genes)
        print(p)
        print(paste0('cluster: ', i)) 
    }
}


options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combinedT, reduction = 'umap',label = T, label.size = 20/.pt, pt.size = 1, group.by = 'Cluster')

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combinedT, reduction = 'umap',label = T, label.size = 20/.pt, pt.size = 1, group.by = 'Loc')

cluster_colors <- c('CD4+TN'='#F57F21','CD4+TCM'='#0FAA07','Tfh_like'='#E41E25','Th2'='#89CC9C',
'ILC2'='#80B4E0','Treg'='#C35338','Th1'='#E05C6E','CD4+CTL1'='#22A3FF',
'CD4+CTL2'='#DCD44E','undefined'='#919191', 'CD4+TMC'='#CE999C')

Idents(combinedT) <- combinedT$seurat_clusters
combinedT <- RenameIdents(combinedT, '2' = 'CD4+TN','0' = 'CD4+TCM', '4' = 'CD4+TMC','8' = 'Tfh_like', '7' = 'Th2', '10' = 'ILC2', 
                         '3' = 'Treg', '6' = 'Th1', '1' = 'CD4+CTL1', '9' = 'CD4+CTL2', 
                         '5' = 'undefined')
options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combinedT, reduction = 'umap',label = T, label.size = 20/.pt, pt.size = 1, repel=T)+
scale_color_manual(values = cluster_colors)



combinedT$model <- combinedT$source
combinedT$model[combinedT$source %in% 'NatMed_asthma'] <- 'Asthma'
combinedT$model[combinedT$source %in% 'GSE179292'] <- 'CRSwNP'

combinedT$new_idents <- Idents(combinedT)

saveRDS(combinedT, file = 'combinedT.rds')