###defined a stem-like Th17 geneset based on public data
library(Seurat)
library(ggplot2)
library(harmony)

load('/data1/02.private/dengyj/analysis/database/mouse_cc_genes_for_cellranger_2.rds')

files <- dir()
files <- files[grep('^GSM', files)]
files <- files[!grepl('doublet', files)]

files

dir()

obj_list <- lapply(files, function(file){
    exp <- Read10X(file)
    doublet_list <- read.table(paste0(file, '_doublets.txt'))
    obj <- CreateSeuratObject(counts = exp[, doublet_list$V1 == 0], min.cells = 3, min.features = 200)
    obj$batch <- file
    obj <- RenameCells(obj, new.names = paste0(colnames(obj), '_', file))
    return(obj)
})


Th17 <- merge(obj_list[[1]], obj_list[2:length(obj_list)])

Th17[['percent.mt']] <- PercentageFeatureSet(Th17, pattern = '^mt-')

options(repr.plot.width=10, repr.plot.height=10)
VlnPlot(Th17, features = 'percent.mt', group.by = 'batch', pt.size = 0)

Th17 <- Th17[, Th17$percent.mt < 10]





Th17 <- NormalizeData(Th17, verbose = T)
Th17 <- FindVariableFeatures(Th17, selection.method = "vst", nfeatures = 2500)
Th17 <- CellCycleScoring(Th17, s.features = s.genes, g2m.features = g2m.genes)
Th17 <- ScaleData(Th17,verbose = FALSE, vars.to.regress = c('S.Score', 'G2M.Score'))
Th17 <- RunPCA(Th17, npcs = 50, verbose = FALSE)


Th17 <- RunHarmony(Th17, 'batch')

options(repr.plot.width=6, repr.plot.height=6)
ElbowPlot(Th17, ndims = 50, reduction = 'harmony')

Th17 <- RunUMAP(Th17, reduction = "harmony", dims = 1:10)

Th17 <- FindNeighbors(Th17, reduction = "harmony", dims = 1:10)

Th17 <- FindClusters(Th17, resolution = 0.2)

options(repr.plot.width=6, repr.plot.height=6)
DefaultAssay(Th17) <- 'RNA'

DimPlot(Th17, reduction = 'umap', label = T)
DimPlot(Th17, reduction = 'umap', group.by= 'Phase')
DimPlot(Th17, reduction = 'umap', group.by= 'batch')


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Th17, reduction = 'umap', cols = c('lightgrey', 'red'), 
            features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'))

###find CD27+ and CD27- Th17
options(repr.plot.width=8, repr.plot.height=8)
FeaturePlot(Th17, reduction = 'umap', cols = c('lightgrey', 'red'), 
            features = c('Cd27'))

options(repr.plot.width=8, repr.plot.height=8)
FeaturePlot(Th17, reduction = 'umap', cols = c('lightgrey', 'red'), 
            features = c('Zbtb16'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Th17, reduction = 'umap', cols = c('lightgrey', 'red'), 
            features = c('Cd3e', 'Cd3d', 'Cd3g', 'Cd4', 'Cd40lg', 'Cd8b1', 'Cd8a'))

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Th17, reduction = 'umap', cols = c('lightgrey', 'red'), 
            features = c('Spi1', 'Tyrobp',  'Ms4a1', 'Gp9', 'Pf4', 'Ptcra', 
                        'Il3ra', 'Cd1c', 'Clec10a', 'Ms4a7'))


markers <- FindMarkers(Th17, pseudocount.use = 0,logfc.threshold = log(1.5), 
                       `ident.1` = c('1', '2'), `ident.2` = c('0', '4', '5'))

saveRDS(markers, file = 'Th17_markers.rds')

GSE121599_Th17_GS <- 
list('CD27P_Th17_up' = 
     rownames(markers[markers$p_val_adj < 0.05 & markers$avg_logFC > log(2), ]), 
    'CD27N_Th17_up' = 
     rownames(markers[markers$p_val_adj < 0.05 & markers$avg_logFC < -log(2), ]))


write_gmt(GSE121599_Th17_GS, filename = 'GSE121599_Th17_genesets.gmt')

saveRDS(Th17, file = 'GSE121599_Th17.rds')
