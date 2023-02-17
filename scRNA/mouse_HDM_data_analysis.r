####to confirm the conservation of Th2 betwenn dust mite and OVA model
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

load('/data1/02.private/dengyj/analysis/database/mouse_cc_genes_for_cellranger_2.rds')

object_list <- lapply(dir('data'), function(x){
    df <-  fread(paste0('data/', x)) 
    df$gene <- make.unique(df$gene)
    df <- data.frame(df, row.names = 'gene')
    sample <- gsub('_rpkms.*', '', x)
    colnames(df) <- paste0(colnames(df), '_', sample)
    object <- CreateSeuratObject(counts = df)
    object$batch <- sample
    object
})

GSE131935 <- merge(object_list[[1]], object_list[2:length(object_list)])

GSE131935[['percent.mt']] <- PercentageFeatureSet(GSE131935, pattern = '^mt-')

GSE131935 <- NormalizeData(GSE131935, verbose = T)
GSE131935 <- FindVariableFeatures(GSE131935,selection.method = "vst", nfeatures = 2000)
GSE131935 <- CellCycleScoring(GSE131935, s.features = s.genes, g2m.features = g2m.genes, 
                                          set.ident = T)
GSE131935$CC.Difference <- GSE131935$S.Score - GSE131935$G2M.Score 

GSE131935 <- ScaleData(GSE131935)#, vars.to.regress = c(names(cc_genes)))
GSE131935 <- RunPCA(GSE131935, npcs = 50)

GSE131935 <- RunHarmony(GSE131935, 'batch')

options(repr.plot.width=7, repr.plot.height=7)
ElbowPlot(GSE131935, ndims = 50, reduction = 'harmony')

GSE131935 <- RunUMAP(GSE131935, reduction = "harmony", dims = 1:10)

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(GSE131935, reduction = 'umap', group.by = 'Phase')

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(GSE131935, reduction = 'umap', group.by = 'batch')

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(GSE131935, features = c('percent.mt', 'nCount_RNA', 'nFeature_RNA'),
            cols = c("lightgrey", 'red'))

GSE131935 <- FindNeighbors(GSE131935, reduction = "harmony", dims = 1:10)

GSE131935 <- FindClusters(GSE131935, resolution = 0.8)

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(GSE131935, reduction = 'umap',label = T, label.size = 20/.pt)

markers_GSE131935 <- FindAllMarkers(GSE131935, pseudocount.use = 0)

openxlsx::write.xlsx(markers_GSE131935, file = 'markers_GSE131935.xlsx')

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(GSE131935, features = c('Cd3e', 'Cd3d', 'Cd3g'),
            cols = c("lightgrey", 'red'))

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(GSE131935, features = c('Foxp3', 'Gata3', 'Cd4', 'Cd8a'),
            cols = c("lightgrey", 'red'))

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(GSE131935, features = c('Il13', 'Il1rl1', 'Gata3', 'Epas1',
                                   'Igfbp7', 'Il5', 'Il4', 'Plac8'),
            cols = c("lightgrey", 'red'))

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(GSE131935, features = c('Tcf7', 'Il7r', 'S1pr1', 'Tmem176b', 'Cd44', 'Anxa1', 'Gpr183', 'Sell', 
                                   'Klf2'),
            cols = c("lightgrey", 'red'))

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(GSE131935, features = c('Cxcr3', 'Ifng', 'Tbx21', 'Gzmk', 'Ccl5', 'Itgb1', 'Itga4'),
            cols = c("lightgrey", 'red'))

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(GSE131935, features = c('Spi1', 'Plek', 'Clec7a','Sirpa'),
            cols = c("lightgrey", 'red'))

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(GSE131935, reduction = 'umap',label = T, label.size = 20/.pt)

Idents(GSE131935) <- GSE131935$seurat_clusters
GSE131935 <- RenameIdents(GSE131935, '3' = 'CD4+TN','2' = 'Th1', '0' = 'Th2','5' = 'IFN_T', '1' = 'Treg', 
                         '4' = 'Treg', 
                         '8' = 'Mki67_T', '7' = 'Myeloid','6' = 'undefined')
options(repr.plot.width=8, repr.plot.height=8)
DimPlot(GSE131935, reduction = 'umap',label = T, label.size = 20/.pt, pt.size = 2)
#guides(color = guide_legend(override.aes = list(size= 7)))

GSE131935$idents <- Idents(GSE131935)

plot_df <- FetchData(GSE131935, vars = c('Il13', 'Il1rl1', 'Gata3', 'Epas1',
                                  'Il5', 'Il4', 'idents'))
plot_df <- plot_df[!plot_df$idents %in% c('Myeloid'), ]
plot_df <- melt(plot_df, id.vars = c('idents'), variable.name = 'features', value.name = 'exp')

plot_df <- ddply(plot_df, c('idents','features'), 
                 summarise, pct = sum(exp>0)/length(exp)*100, exp = mean(exp))

plot_df <- ddply(plot_df, c('features'), transform, exp = scale(exp))

options(repr.plot.width=6, repr.plot.height=4)
ggplot(plot_df, aes(features, idents, fill=exp, size= pct))+
geom_point(shape=21)+
scale_fill_gradientn(colours = c('blue', 'lightgrey', 'red'))




OVA_Th2 <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/basic_analysis/Hif2a_Th2_2.rds')

OVA_Th2<- OVA_Th2[rowSums(OVA_Th2$RNA@data)>0, ]

OVA_Th2$model <- 'OVA'

HDM_Th2 <- GSE131935[, Idents(GSE131935) %in% c("Th2")]

HDM_Th2<- HDM_Th2[rowSums(HDM_Th2$RNA@data)>0, ]

HDM_Th2$model <- 'HDM'

common_features <- intersect(rownames(OVA_Th2), rownames(HDM_Th2))

combined_Th2 <- merge(OVA_Th2[common_features, ], HDM_Th2[common_features, ])

table(combined_Th2$batch)

DefaultAssay(combined_Th2) <- 'RNA'
combined_Th2_list <- SplitObject(combined_Th2, split.by = 'batch')
combined_Th2_list <- lapply(combined_Th2_list, function(obj){
    obj <- FindVariableFeatures(obj)
    obj
})

anchors_features <- SelectIntegrationFeatures(combined_Th2_list, nfeatures = 2000)

anchors <- FindIntegrationAnchors(object.list = combined_Th2_list, dims = 1:20, 
                                  anchor.features = anchors_features, k.filter = 40)

combined_Th2 <- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(combined_Th2) <- 'integrated'
combined_Th2 <- ScaleData(combined_Th2)#, vars.to.regress = names(Go_list))
combined_Th2 <- RunPCA(combined_Th2, npcs = 50, verbose = T)

options(repr.plot.width=6, repr.plot.height=6)
ElbowPlot(combined_Th2, ndims = 50)

combined_Th2 <- RunTSNE(combined_Th2, reduction = "pca", dims = setdiff(1:20, c()), verbose = F)

#combined_Th2 <- RunUMAP(combined_Th2, reduction = "pca", dims = setdiff(1:15, c()), verbose = F)

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(combined_Th2[, combined_Th2$model %in% c('OVA')], reduction = 'tsne', label = T, group.by = 'seurat_clusters')

pca_df <- Embeddings(combined_Th2, reduction = 'pca')


OVA_pca_df <- data.frame(pca_df[colnames(combined_Th2)[combined_Th2$model %in%  c('OVA')], ])
OVA_pca_df <- OVA_pca_df[, paste0('PC_', 1:15)]
OVA_pca_df$orig_idents <- factor(combined_Th2$seurat_clusters[combined_Th2$model %in%  c('OVA')])

train_cells <- lapply(unique(OVA_pca_df$orig_idents), function(id){
    df_iter <- OVA_pca_df[OVA_pca_df$orig_idents %in% id, ]
    cells <- rownames(df_iter)
    set.seed(1234)
    cells_selected <- sample(cells, size = floor(length(cells)*0.7))
    cells_selected
})
train_cells <- unlist(train_cells)

set.seed(123)
fit.forest <- randomForest(orig_idents~., data=OVA_pca_df[train_cells, ],  importance=TRUE)
forest.pred <- predict(fit.forest, OVA_pca_df[setdiff(rownames(OVA_pca_df), train_cells), ])
forest.perf <- table(OVA_pca_df[setdiff(rownames(OVA_pca_df), train_cells), 'orig_idents'],
                     forest.pred,dnn=c("Actual", "Predicted"))

roc_result <- roc(as.numeric(OVA_pca_df[setdiff(rownames(OVA_pca_df), train_cells), 'orig_idents']), 
                  as.numeric(forest.pred))

roc_result

HDM_pca_df <- data.frame(pca_df[colnames(combined_Th2)[combined_Th2$model %in%  c('HDM')], ])

HDM_forest.pred <- data.frame(HDM_rf_pred = predict(fit.forest, HDM_pca_df))


combined_Th2 <- AddMetaData(combined_Th2, metadata = HDM_forest.pred)
combined_Th2$HDM_rf_pred <- as.character(combined_Th2$HDM_rf_pred)

combined_Th2$combined_ID <- combined_Th2$seurat_clusters
combined_Th2$combined_ID[combined_Th2$model %in% 'HDM'] <- combined_Th2$HDM_rf_pred[combined_Th2$model %in% 'HDM']

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

plot_df <- cbind(FetchData(combined_Th2, vars = c('combined_ID', 'model')), Embeddings(combined_Th2, reduction = 'tsne'))
options(repr.plot.width=12, repr.plot.height=6)
ggplot(plot_df, aes(tSNE_1, tSNE_2, color = combined_ID))+
geom_point()+
# stat_ellipse(data = plot_df[plot_df$model %in% c('mouse'), ], mapping = aes(group = combined_ID),
#              level = 0.7, linetype='dashed', size = 1, show.legend =F, color = 'red')+
#scale_color_manual(values = cluster_colors)+
facet_wrap(.~model)
Idents(combined_Th2) <- combined_Th2$combined_ID
table(combined_Th2$combined_ID, combined_Th2$model)

####
combined_Th2$final_ID <- as.character(combined_Th2$combined_ID)
combined_Th2$final_ID[combined_Th2$combined_ID %in% c('3')] <- '0'###stem-like Th2
combined_Th2$final_ID[combined_Th2$combined_ID %in% c('0')] <- '1'#pathogenic Th2
combined_Th2$final_ID[combined_Th2$combined_ID %in% c('1')] <- '3'
#DefaultAssay(combined_Th2) <- 'symbol'

Idents(combined_Th2) <- combined_Th2$final_ID

plot_df <- cbind(FetchData(combined_Th2, vars = c('final_ID', 'model')), Embeddings(combined_Th2, reduction = 'tsne'))
options(repr.plot.width=12, repr.plot.height=6)
ggplot(plot_df, aes(tSNE_1, tSNE_2, color = final_ID))+
geom_point()+

scale_color_manual(values = cluster_colors)+
facet_wrap(.~model)



combined_Th2_conserved_markers_list <- lapply(unique(combined_Th2$final_ID), function(id){
    if(id != '5'){
        marker <- FindConservedMarkers(combined_Th2, ident.1 = id, 
                                       grouping.var = 'model',assay = 'RNA',# group.by='final_ID', 
                                       pseudocount.use=0)   
        marker$cluster <- id
        marker        
    }
})

names(combined_Th2_conserved_markers_list) <- unique(combined_Th2$combined_ID)
saveRDS(combined_Th2_conserved_markers_list, file = 'combined_Th2_conserved_markers_list.rds')

wb <- createWorkbook()
for(i in names(combined_Th2_conserved_markers_list)){
    addWorksheet(wb, i)
    writeData(wb, sheet = i, x = combined_Th2_conserved_markers_list[[i]], rowNames = T)
}
saveWorkbook(wb = wb, file ='combined_Th2_conserved_markers_list.xlsx',overwrite = T)

plot_df <- FetchData(combined_Th2, 
                     vars = c('Slamf6','Rel','Cd200','Egr2',
                              'Il1rl1','Gata3','Epas1','Ier3','Cd27','Ctla4','Gimap4',
                              'Cxcr3', 'final_ID', 'model'))
###'2' represents Treg-like Th2                              
plot_df <- plot_df[plot_df$final_ID %in% c('0','1','2'), ]
plot_df <- melt(plot_df, id.vars = c('final_ID', 'model'), variable.name = 'features', value.name = 'exp')


plot_df <- ddply(plot_df, c('final_ID', 'model','features'), 
                 summarise, pct = sum(exp>0)/length(exp)*100, exp = mean(exp))

plot_df <- ddply(plot_df, c('model','features'), transform, exp = scale(exp))


options(repr.plot.width=8, repr.plot.height=4)
ggplot(plot_df, aes(features, final_ID, fill=exp, size= pct))+
geom_point(shape=21)+
facet_wrap(.~model, scale='free_x')+
scale_fill_gradientn(colours = c('blue', 'lightgrey', 'red'))
saveRDS(GSE131935, file = 'GSE131935.rds')

saveRDS(combined_Th2, file = 'combined_Th2.rds')
