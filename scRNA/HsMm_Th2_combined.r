####to confirm the conservation of Th2 betwenn human and mouse
library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)
library(openxlsx)
library(future)
library(org.Mm.eg.db)
library(SCopeLoomR)
library(plyr)

Mm_Th2 <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/basic_analysis/Hif2a_Th2_2.rds')

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(Mm_Th2, features = c('Lef1','Tcf7', 'Sell', 'Ccr7'),
            cols = c("lightgrey", 'red'))

Hs_Th2 <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/public_data/Hs_combined/combined_Th2.rds')





genes_tsv <- read.table('/data1/01.project/Hif2a//ABFC20200274-060709（ctm_hif2a）//RNA/ABFC20200274_10X-scRNA_result/2.Cellranger/W1-10X5/filtered_feature_bc_matrix/features.tsv.gz')
colnames(genes_tsv) <- c('ensembl', 'symbol')

genes_tsv$symbol <- make.unique(genes_tsv$symbol)

Hs_to_Mm <- 
data.table::fread('/data1/02.private/dengyj/analysis/database/cell_cycle/Whitfield/mart_export_Homologues_104.txt')
Hs_to_Mm <- data.frame(Hs_to_Mm)

head(Hs_to_Mm)

Hs_to_Mm <- Hs_to_Mm[!is.na(Hs_to_Mm$`Mouse.orthology.confidence..0.low..1.high.`), ]
Hs_to_Mm <- Hs_to_Mm[Hs_to_Mm$`Mouse.orthology.confidence..0.low..1.high.` == 1, ]

exp <- Hs_Th2$RNA@data

exp_filtered <- exp[rownames(exp) %in% Hs_to_Mm$Gene.name, ]
#exp_filtered <- exp_filtered[rowSums(exp_filtered) > 0, ]

new_names <- lapply(rownames(exp_filtered), function(x){
    Mm_ensembl <- unique(Hs_to_Mm[Hs_to_Mm$Gene.name %in% x, 'Mouse.gene.stable.ID'])
    if(length(Mm_ensembl) > 0){
        return(Mm_ensembl)
    }else{
        return(NULL)
    }
})

exp_filtered_Mm <- exp_filtered[sapply(new_names, length) == 1, ]

Hs_mapping_Mm <- data.frame(Hs_symbol = rownames(exp_filtered_Mm), 
                            Mm_ensembl = unlist(new_names[sapply(new_names, length) ==1]))

rownames(exp_filtered_Mm) <- unlist(new_names[sapply(new_names, length) ==1])

sum(duplicated(rownames(exp_filtered_Mm)))

exp_filtered_Mm <- exp_filtered_Mm[!rownames(exp_filtered_Mm) %in%#####去除重复的基因
                                   rownames(exp_filtered_Mm)[duplicated(rownames(exp_filtered_Mm))], ]

Hs_mapping_Mm <- Hs_mapping_Mm[!duplicated(Hs_mapping_Mm$Mm_ensembl), ]

sum(duplicated(rownames(exp_filtered_Mm)))

exp_filtered_Mm[1:10, 1:10]

get_id_changed <- function(
  orig_id, 
  mapping, 
  orig_var = 'V2', 
  new_var = 'V1'){
  result <- sapply(orig_id, function(i){
    if(i %in% mapping[, orig_var]){
      val <- mapping[mapping[, orig_var] %in% i, new_var][1]
    }else{
      loc <- stringr::str_extract_all(i, '\\.\\d*$')
      loc <- as.numeric(gsub('\\.', '', loc))
      i_cp <- gsub('\\.\\d*$','',  i)
      val <- mapping[mapping[, orig_var] %in% i_cp, new_var][(loc+1)]
    }
    return(val)}
  )
  return(result)
}

Mm_Th2_ensembl_exp <- Mm_Th2$RNA@data

Mm_Th2_ensembl_id <- 
get_id_changed(rownames(Mm_Th2$RNA@data), mapping = genes_tsv, orig_var = 'symbol', new_var = 'ensembl')

rownames(Mm_Th2_ensembl_exp) <- Mm_Th2_ensembl_id

Ms_symbol_ensembl <- data.frame(Ms_symbol = rownames(Mm_Th2$RNA@data), Mm_ensembl = 
                               Mm_Th2_ensembl_id)

sum(duplicated(rownames(Mm_Th2_ensembl_exp)))

Mm_Th2_ensembl_exp[1:10,1:10]



common_features <- intersect(rownames(exp_filtered_Mm), rownames(Mm_Th2_ensembl_exp))

length(common_features)

Hs_mapping_Mm <- Hs_mapping_Mm[Hs_mapping_Mm$Mm_ensembl %in% common_features, ]
Ms_symbol_ensembl <- Ms_symbol_ensembl[Ms_symbol_ensembl$Mm_ensembl %in% common_features, ]

mapping_df <- merge(Hs_mapping_Mm, Ms_symbol_ensembl)

head(mapping_df)



new_Hs_Th2 <- CreateSeuratObject(counts = exp_filtered_Mm[common_features, ], meta.data = Hs_Th2@meta.data)
new_Hs_Th2$species <- 'human'
new_Hs_Th2$batch_to_reduce <- new_Hs_Th2$source
new_Hs_Th2$orig_idents <- as.character(Idents(Hs_Th2))

Mm_Th2_ensembl <- CreateSeuratObject(counts = Mm_Th2_ensembl_exp[common_features, ], meta.data = Mm_Th2@meta.data)
Mm_Th2_ensembl$species <- 'mouse'
Mm_Th2_ensembl$batch_to_reduce <- Mm_Th2_ensembl$batch
Mm_Th2_ensembl$orig_idents <- as.character(Mm_Th2_ensembl$seurat_clusters)

dim(Mm_Th2_ensembl)

HsMm_Th2 <- merge(Mm_Th2_ensembl, new_Hs_Th2)

HsMm_Th2_symbol_mat <- HsMm_Th2$RNA@data

rownames(HsMm_Th2_symbol_mat) <- mapping_df[match(rownames(HsMm_Th2_symbol_mat), mapping_df$Mm_ensembl), 'Ms_symbol']

HsMm_Th2[['symbol']] <- CreateAssayObject(counts = HsMm_Th2_symbol_mat)



DefaultAssay(HsMm_Th2) <- 'RNA'
HsMm_Th2_list <- SplitObject(HsMm_Th2, split.by = 'batch_to_reduce')
HsMm_Th2_list <- lapply(HsMm_Th2_list, function(obj){
    obj <- FindVariableFeatures(obj)
    obj
})

anchors_features <- SelectIntegrationFeatures(HsMm_Th2_list, nfeatures = 2000)

anchors <- FindIntegrationAnchors(object.list = HsMm_Th2_list, dims = 1:20, 
                                  anchor.features = anchors_features, k.filter = 50)
HsMm_Th2 <- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(HsMm_Th2) <- 'integrated'
HsMm_Th2 <- ScaleData(HsMm_Th2)#, vars.to.regress = names(Go_list))
HsMm_Th2 <- RunPCA(HsMm_Th2, npcs = 50, verbose = T)

options(repr.plot.width=6, repr.plot.height=6)
ElbowPlot(HsMm_Th2, ndims = 50)

HsMm_Th2 <- RunTSNE(HsMm_Th2, reduction = "pca", dims = setdiff(1:15, c()), verbose = F)

HsMm_Th2 <- FindNeighbors(HsMm_Th2, reduction = "pca", dims = setdiff(1:15, c()))

HsMm_Th2 <- FindClusters(HsMm_Th2, resolution = 0.5)

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(HsMm_Th2, reduction = "tsne", label = T, group.by = 'seurat_clusters')

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(HsMm_Th2[, HsMm_Th2$species %in% c('mouse')], reduction = "tsne", label = T, group.by = 'orig_idents')

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(HsMm_Th2, reduction = "tsne", features = c('Cxcr6', 'Cxcr5','Bcl6', 'Pdcd1'))

colors <- c('#F57F21', '#919191','#0FAA07', '#C35338', '#89CC9C', '#80B4E0',
 '#DCD44E', '#E05C6E', '#64446C', '#22A3FF', '#E41E25','#CE999C',  
 '#136433', '#9C56A3', '#FE3DF6', '#9F8CFF')



library(randomForest) 

library(pROC)


pca_df <- Embeddings(HsMm_Th2, reduction = 'pca')

mouse_pca_df <- data.frame(pca_df[colnames(HsMm_Th2)[HsMm_Th2$species %in%  c('mouse')], ])

mouse_pca_df <- mouse_pca_df[, paste0('PC_', 1:15)]
mouse_pca_df$orig_idents <- factor(HsMm_Th2$orig_idents[HsMm_Th2$species %in%  c('mouse')])

train_cells <- lapply(unique(mouse_pca_df$orig_idents), function(id){
    df_iter <- mouse_pca_df[mouse_pca_df$orig_idents %in% id, ]
    cells <- rownames(df_iter)
    set.seed(1234)
    cells_selected <- sample(cells, size = floor(length(cells)*0.7))
    cells_selected
})
train_cells <- unlist(train_cells)

set.seed(123)
fit.forest <- randomForest(orig_idents~., data=mouse_pca_df[train_cells, ],  importance=TRUE)

forest.pred <- predict(fit.forest, mouse_pca_df[setdiff(rownames(mouse_pca_df), train_cells), ])

forest.perf <- table(mouse_pca_df[setdiff(rownames(mouse_pca_df), train_cells), 'orig_idents'],
                     forest.pred,dnn=c("Actual", "Predicted"))

forest.perf

roc_result <- roc(as.numeric(mouse_pca_df[setdiff(rownames(mouse_pca_df), train_cells), 'orig_idents']), 
                  as.numeric(forest.pred))

roc_result




human_pca_df <- data.frame(pca_df[colnames(HsMm_Th2)[HsMm_Th2$species %in%  c('human')], ])

Hs_forest.pred <- data.frame(Hs_rf_pred = predict(fit.forest, human_pca_df))

HsMm_Th2 <- AddMetaData(HsMm_Th2, metadata = Hs_forest.pred)
HsMm_Th2$Hs_rf_pred <- as.character(HsMm_Th2$Hs_rf_pred)

HsMm_Th2$combined_ID <- HsMm_Th2$orig_idents
HsMm_Th2$combined_ID[HsMm_Th2$species %in% 'human'] <- HsMm_Th2$Hs_rf_pred[HsMm_Th2$species %in% 'human']

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

plot_df <- cbind(FetchData(HsMm_Th2, vars = c('combined_ID', 'species')), Embeddings(HsMm_Th2, reduction = 'tsne'))
options(repr.plot.width=12, repr.plot.height=6)
ggplot(plot_df, aes(tSNE_1, tSNE_2, color = combined_ID))+
geom_point()+
# stat_ellipse(data = plot_df[plot_df$species %in% c('mouse'), ], mapping = aes(group = combined_ID),
#              level = 0.7, linetype='dashed', size = 1, show.legend =F, color = 'red')+
scale_color_manual(values = cluster_colors)+
facet_wrap(.~species)+

#ggsave('/data1/02.private/dengyj/analysis/Hif2a/plot/HsMm_Th2_tsne.svg' ,width = 12, height = 6)

HsMm_Th2$final_ID <- as.character(HsMm_Th2$combined_ID)
HsMm_Th2$final_ID[HsMm_Th2$combined_ID %in% c('3')] <- '0'###stem-like Th2
HsMm_Th2$final_ID[HsMm_Th2$combined_ID %in% c('0')] <- '1'#pathogenic Th2
HsMm_Th2$final_ID[HsMm_Th2$combined_ID %in% c('1')] <- '3'
#DefaultAssay(HsMm_Th2) <- 'symbol'

Idents(HsMm_Th2) <- HsMm_Th2$final_ID

saveRDS(HsMm_Th2, file = 'HsMm_Th2.rds')

plot_df <- cbind(FetchData(HsMm_Th2, vars = c('final_ID', 'species')), Embeddings(HsMm_Th2, reduction = 'tsne'))
options(repr.plot.width=12, repr.plot.height=6)
ggplot(plot_df, aes(tSNE_1, tSNE_2, color = final_ID))+
geom_point()+
# stat_ellipse(data = plot_df[plot_df$species %in% c('mouse'), ], mapping = aes(group = final_ID),
#              level = 0.7, linetype='dashed', size = 1, show.legend =F, color = 'red')+
scale_color_manual(values = cluster_colors)+
facet_wrap(.~species)

HsMm_Th2_conserved_markers_list <- lapply(unique(HsMm_Th2$final_ID), function(id){
    marker <- FindConservedMarkers(HsMm_Th2, ident.1 = id, 
                                   grouping.var = 'species',assay = 'symbol',# group.by='final_ID', 
                                   pseudocount.use=0)   
    marker$cluster <- id
    marker
})


names(HsMm_Th2_conserved_markers_list) <- unique(HsMm_Th2$final_ID)
saveRDS(HsMm_Th2_conserved_markers_list, file = 'HsMm_Th2_conserved_markers_list.rds')



wb <- createWorkbook()
for(i in names(HsMm_Th2_conserved_markers_list)){
    addWorksheet(wb, i)
    writeData(wb, sheet = i, x = HsMm_Th2_conserved_markers_list[[i]], rowNames = T)
}
saveWorkbook(wb = wb, file ='HsMm_Th2_conserved_markers_list.xlsx',overwrite = T)



DefaultAssay(HsMm_Th2) <- 'symbol'
plot_df <- FetchData(HsMm_Th2, 
                     vars = c('Egr2','Tcf7','Il21','Cd200','Il1rl1','Gata3','Epas1','Pparg',
                              'Cd27','Lag3','Ctla4','Tbc1d4',
                               'final_ID', 'species'))
plot_df <- plot_df[plot_df$final_ID %in% c('0','1','2'), ]
plot_df <- melt(plot_df, id.vars = c('final_ID', 'species'), variable.name = 'features', value.name = 'exp')


plot_df <- ddply(plot_df, c('final_ID', 'species','features'), 
                 summarise, pct = sum(exp>0)/length(exp)*100, exp = mean(exp))


plot_df <- ddply(plot_df, c('species','features'), transform, exp = scale(exp))
Mm_lvls <- levels(plot_df$features)
Hs_lvls <- mapping_df[match(Mm_lvls, mapping_df$Ms_symbol), 'Hs_symbol']
combined_lvls <- unlist(lapply(1:length(Mm_lvls), function(iter){
    c(Hs_lvls[iter],Mm_lvls[iter])
}))
plot_df$features <- as.character(plot_df$features)

features_selected <- plot_df[plot_df$species %in% c('human'), 'features']
plot_df[plot_df$species %in% c('human'), 'features'] <- 
mapping_df[match(features_selected, mapping_df$Ms_symbol), 'Hs_symbol']
plot_df$features <- factor(plot_df$features, levels = combined_lvls)

options(repr.plot.width=8, repr.plot.height=4)
ggplot(plot_df, aes(features, final_ID, fill=exp, size= pct))+
geom_point(shape=21)+
facet_wrap(.~species, scale='free_x')+
scale_fill_gradientn(colours = c('blue', 'lightgrey', 'red'))+
theme(axis.title = element_blank(), 
      strip.text.x=element_text(size = 12.5),
      legend.position = 'right',
      legend.text = element_text(size = 12.5, family = 'sans'), 
      legend.title = element_text(size = 12.5, family = 'sans'), 
      plot.title = element_text(size = 12.5, face = 'italic',hjust=0.5), 
      axis.text.y = element_text(hjust = 1, 
                                 size = 12.5, family = 'sans', face = 'plain'),
     axis.text.x = element_text(hjust =1, angle =90,size = 12.5, family = 'sans', face = 'italic'), 
      axis.line = element_blank(), 
     panel.background = element_rect(color = NA, fill = 'white'), 
     panel.border = element_rect(fill = NA, colour = 'black'))+
guides(size = guide_legend(title = "% Cells",override.aes = list(fill = 'black')), 
       fill = guide_colorbar(order = 1,title = "scaled Exp"))
