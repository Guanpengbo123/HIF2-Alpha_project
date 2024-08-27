library(Seurat)
library(ggplot2)
library(future)
library(openxlsx)
library(harmony)

library(Matrix)

library(DoubletFinder)

library(rlang)

library(plyr)

library(patchwork)

library(reshape2)

library(ggforce)
library(VISION)

library(ComplexHeatmap)

source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')

color_used <- c('#E5D2DD', '#53A85F', '#F1BB72',  '#D6E7A3', '#57C3F3', '#476D87',
'#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
'#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
'#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
'#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#F3B1A0',
'#968175','#F57F21', '#919191','#0FAA07', '#C35338', '#89CC9C', '#80B4E0',
 '#DCD44E', '#E05C6E', '#64446C', '#22A3FF', '#E41E25','#CE999C',  
 '#136433', '#9C56A3', '#FE3DF6', '#9F8CFF'
)


source('/data1/02.private/dengyj/analysis/Hif2a/plot/support/plot_helper.r')

get_sem <- function(x){
    sem <- sd(x)/sqrt(length(x))
    val <- c((mean(x) -sem), (mean(x)+sem))
    names(val) <- c('ymin', 'ymax')
    return(val)
}

#GSE175930_Tcell <- readRDS('../Hs_EoE_SciImm_2021/GSE175930_Tcell.rds')

GSE179292_2 <- readRDS('../Hs_nasal_polyp_SciImm_2021/GSE179292_2.rds')

NatMed_CD4_T <- readRDS('../Hs_asthma_NatMed_2019/CD4_T.rds')

CD4T_ILC2 <- readRDS('../Hs_CRS_NI_2022/CD4T_ILC2.rds')

SI_CD4_T <- readRDS('../Hs_asthma_SI_2023/GSE193816_AA_CD4T_3.rds')

#GSE175930_Tcell$source <- 'GSE175930'

GSE179292_2$source <- 'GSE179292'

NatMed_CD4_T$source <- 'NatMed_asthma'
NatMed_CD4_T$batch <- NatMed_CD4_T$Loc

CD4T_ILC2$source <- 'NI_CRS'

SI_CD4_T$source <- 'SI'
SI_CD4_T$batch <- paste0(SI_CD4_T$id,'_', SI_CD4_T$sample)

common_features <- intersect(rownames(NatMed_CD4_T), 
                             union(rownames(GSE179292_2), union(rownames(CD4T_ILC2), rownames(SI_CD4_T))))
# common_features <- intersect(common_features, rownames(CD4T_ILC2))
# common_features <- intersect(common_features, rownames(SI_CD4_T))

combinedT <- merge(GSE179292_2[, !GSE179292_2$Identity %in% c('non_Lym')], 
                  #GSE175930_Tcell[common_features, GSE175930_Tcell$Identity %in% c('CD4T')], 
                      list(NatMed_CD4_T, CD4T_ILC2, SI_CD4_T))

combinedT[['common']] <- CreateAssayObject(counts = combinedT$RNA@counts[common_features, ])

combinedT$common@data <- combinedT$RNA@data[rownames(combinedT$common), ]

orig_idents <- readRDS('orig_idents.rds')

combinedT$model <- combinedT$source
combinedT$model[combinedT$source %in% 'NI_CRS'] <- 'CRSwNP'
combinedT$model[combinedT$source %in% 'NatMed_asthma'] <- 'Asthma'
combinedT$model[combinedT$source %in% 'SI'] <- 'Asthma'
combinedT$model[combinedT$source %in% 'GSE179292'] <- 'CRSwNP'


combinedT <- AddMetaData(combinedT, orig_idents, col.name = 'orig_idents')



options(future.globals.maxSize = 10 * 1024^3)
#plan('multisession', workers = 10)
DefaultAssay(combinedT) <- 'common'
combinedT <- FindVariableFeatures(combinedT)
combinedT <- ScaleData(combinedT)#, vars.to.regress = c('nCount_RNA'))
#plan('sequential')


combinedT <- RunPCA(combinedT, npcs = 50, verbose = F)

unique(combinedT$source)

combinedT <- RunHarmony(combinedT, c('source', 'batch'),assay.use = 'common')

options(repr.plot.width=6, repr.plot.height=6)
ElbowPlot(combinedT, ndims = 50, reduction = 'harmony')

DefaultAssay(combinedT) <- 'common'

combinedT <- RunUMAP(combinedT, reduction = "harmony", dims = 1:30, verbose = F)

#combinedT <- RunTSNE(combinedT, reduction = "pca", dims = 1:10)

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combinedT, reduction = 'umap', group.by = 'Phase')
options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combinedT, reduction = 'umap', group.by = 'source')


options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combinedT, reduction = 'umap', group.by = 'source', sizes.highlight = 0.1)


Th2_cells <- readRDS('../Hs_asthma_SI_2023/Th2_cells.rds')

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combinedT, reduction = 'umap', cells.highlight = Th2_cells, sizes.highlight = 0.1)


options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combinedT, reduction = 'umap', group.by = 'sample')+
scale_color_manual(values = color_used)

options(repr.plot.width=8, repr.plot.height=8)
DimPlot(combinedT, reduction = 'umap', group.by = 'orig_idents')+
scale_color_manual(values = color_used)



DefaultAssay(combinedT) <- 'RNA'


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('IL13', 'IL1RL1', 'GATA3', 'EPAS1',
                                   'CD200R1', 'IL5', 'IL4', 'PTGDR2', 'IL10'),
            cols = c("lightgrey", 'red'))                        

####Th2， ILC2
options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('TCN1','HPGDS','CD109','PDCD1', 'CD200', 'TOX2', 'FOXP3', 'XCL1', 'CD3E'),
            cols = c("lightgrey", 'red'))   


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, reduction = 'umap',features = c('TCF7', 'LEF1', 'SELL', 'BATF'),
            cols = c("lightgrey", 'red'))        


options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('KLRB1'),
            cols = c("lightgrey", 'red'))        



DefaultAssay(combinedT) <- 'RNA'

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('CD4', 'CD40LG', 'CD8A', 'CD8B'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('IL17A','TCF7','GZMB','FOXP3'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('MKI67'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('CD3D', 'CD3E', 'CD3G'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('S1PR1', 'KLF2', 'NR4A2'),ncol = 2,
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('SLC4A10', 'RORC', 'NCR3', 'ZBTB16'),
            cols = c("lightgrey", 'red'))              


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('TBX21', 'IFNG', 'CCR6', 'RORC'),
            cols = c("lightgrey", 'red'))        


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('IL17A', 'IL17F', 'CCR6', 'RORC'),
            cols = c("lightgrey", 'red'))        


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('TMEM176A', 'TMEM176B', 'IL17RE','IL23R'),
            cols = c("lightgrey", 'red'))        


options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combinedT, features = c('FOXP3', 'CTLA4', 'IL2RA'),
            cols = c("lightgrey", 'red'))        


options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('PTPRC'),
            cols = c("lightgrey", 'red')) 

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, reduction = 'umap',features = c('CD44', 'FAS', 'ANXA2', 'PASK', 
                                         'IL7R', 'GPR183', 'AHNAK', 'ANXA1'),#ncol=4,
            cols = c("lightgrey", 'red')) 

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, reduction = 'umap',features = c('PDCD1', 'BCL6', 'TOX2', 'IL6ST'),
            cols = c("lightgrey", 'red')) 

combinedT=readRDS('combinedT.rds')

dim(combinedT)

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combinedT, reduction = 'umap', cols = c('lightgrey', 'red'),
                 #min.cutoff = 'q9',
                 features = c('CD27','CTLA4','TIGIT', 'LAG3'
                             ))#Eosinophil

options(repr.plot.width=16, repr.plot.height=16)
FeaturePlot(combinedT, reduction = 'umap', cols = c('lightgrey', 'red'),
                 #min.cutoff = 'q9',
                 features = c('FOXP1','FOXP3','CCR5', 'IKZF2'
                             ))

options(repr.plot.width=16, repr.plot.height=8)
FeaturePlot(combinedT, reduction = 'umap', cols = c('lightgrey', 'red'),
                 #min.cutoff = 'q9',
                 features = c('TNFRSF9','TNFRSF18'
                             ))

options(repr.plot.width=8, repr.plot.height=8)
FeaturePlot(combinedT, reduction = 'umap', cols = c('lightgrey', 'red'),
                 #min.cutoff = 'q9',
                 features = c('FOXP3'
                             ))


options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('CXCR5'),
            cols = c("lightgrey", 'red'))        


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('PLA2G16', 'PTGS2','HPGDS'),
            cols = c("lightgrey", 'red'))        


options(repr.plot.width=16, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('KLRB1','PTGDR2'),
            cols = c("lightgrey", 'red'))        


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('TCF7','SELL', 'S1PR1', 'CCR7'),
            cols = c("lightgrey", 'red'))        


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('GZMB', 'GZMK', 'PRF1', 'CX3CR1'),
            cols = c("lightgrey", 'red'))        

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, features = c('CD8A', 'CD8B', 'SLC4A10', 'NCR3'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT[, !colnames(combinedT) %in% colnames(SI_CD4_T)],
            features = c('CD8A', 'CD8B', 'SLC4A10', 'NCR3'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=16, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('PTGDR2', 'KLRB1'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('MS4A1'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('ACTN1'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('CD109'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('FAS'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('ANXA2'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('AHNAK'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('EPAS1'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('CXCR6'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('CD27'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('PTGDR2'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('KLRG1'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('IL17RB'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('IL1RL1'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('S1PR1'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('KLRB1'),
            cols = c("lightgrey", 'red'))  

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('SATB1'),
            cols = c("lightgrey", 'red'))  



DefaultAssay(combinedT) <- 'common'
combinedT <- FindNeighbors(combinedT, reduction = "harmony", dims = 1:30)
DefaultAssay(combinedT) <- 'RNA'

DefaultAssay(combinedT) <- 'common'
combinedT <- FindClusters(combinedT, resolution = seq(0.1,1.2, 0.05), verbose = F)
DefaultAssay(combinedT) <- 'RNA'

options(repr.plot.width=6, repr.plot.height=6)
for(i in seq(0.1,1.2, 0.05)){
    res <- paste0('common_snn_res.', i)
    p <- DimPlot(combinedT, label = T, group.by= res)+
    labs(title = res)
    print(p)
}

options(future.globals.maxSize = 10 * 1024^3)
DefaultAssay(combinedT) <- 'RNA'
Idents(combinedT) <- combinedT$common_snn_res.0.6
plan('multisession', workers = 10)
markers_combinedT <- FindAllMarkers(combinedT, logfc.threshold = log(1.5), pseudocount.use = 0)
plan('sequential')


write.xlsx(markers_combinedT, file = 'markers_combinedT.xlsx')

combinedT <- ScaleData(combinedT, features = rownames(combinedT))

p <- DoHeatmap(combinedT, size = 2,label=T, #group.colors = cluster_colors,
               features = unique(markers_combinedT[markers_combinedT$p_val_adj <0.05 &
                                               markers_combinedT$avg_logFC > log(2) &
                                               markers_combinedT$pct.1 > 0.3, 'gene']) ) + 
  theme(legend.text = element_text(size=6), 
        legend.title = element_text(size = 6), 
        axis.text.y = element_blank())
options(repr.plot.width=8, repr.plot.height=8)
p

options(future.globals.maxSize = 10 * 1024^3)
plan('multisession', workers = 10)
markers <- FindMarkers(combinedT, pseudocount.use = 0, `ident.1` = '9',`ident.2` = '8', 
                       group.by = 'common_snn_res.0.65')
plan('sequential')
head(markers[markers$p_val_adj < 0.05 & markers$avg_logFC > 0,], 50)

head(markers[markers$p_val_adj < 0.05 & markers$avg_logFC < 0,], 50)

options(repr.plot.width=6, repr.plot.height=6)
DimPlot(combinedT, label = T, cells.highlight = colnames(combinedT)[combinedT$`common_snn_res.0.6`=='12'])


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT, reduction = 'umap',features = c('TCF7', 'MAL', 'SELL', 'KLF2'),
            cols = c("lightgrey", 'red'))        


options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, reduction = 'umap',features = c('S1PR1'),
            cols = c("lightgrey", 'red'))        

options(repr.plot.width=8, repr.plot.height=8)#
FeaturePlot(combinedT, features = c('AREG'),
            cols = c("lightgrey", 'red'))

DefaultAssay(combinedT) <- 'common'
options(repr.plot.width=6, repr.plot.height=4)#
DotPlot(combinedT, #reduction = 'umap',
        features = c('HIF1A', 'EPAS1'),
            cols = c("lightgrey", 'red'))                        




combinedT_S1 <- combinedT[, combinedT$common_snn_res.0.6 == '6']
DefaultAssay(combinedT_S1) <- 'common'

combinedT_S1 <- FindNeighbors(combinedT_S1, reduction = "harmony", dims = 1:15)

combinedT_S1 <- FindClusters(combinedT_S1, resolution = 0.2, verbose = F)

options(repr.plot.width=6, repr.plot.height=6)
 DimPlot(combinedT_S1, label = T)

combinedT_S1$Identity <- as.character(combinedT_S1$seurat_clusters)
combinedT_S1$Identity[combinedT_S1$Identity %in% c('1')]<- 'Th2'
combinedT_S1$Identity[combinedT_S1$Identity %in% c('0')]<- 'ILC2'
Idents(combinedT_S1) <- combinedT_S1$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(combinedT_S1, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)

combinedT_S2 <- combinedT[, combinedT$common_snn_res.0.6 == '7']
DefaultAssay(combinedT_S2) <- 'common'

combinedT_S2 <- FindNeighbors(combinedT_S2, reduction = "harmony", dims = 1:15)

combinedT_S2 <- FindClusters(combinedT_S2, resolution = 0.5, verbose = F)

options(repr.plot.width=6, repr.plot.height=6)
 DimPlot(combinedT_S2, label = T)


options(repr.plot.width=16, repr.plot.height=16)#
FeaturePlot(combinedT_S2, reduction = 'umap',features = c('CXCR5', 'PDCD1', 'TOX2', 'IL4', 'IL5', 'GATA3',
                                                         'IL13', 'IL1RL1'),
            cols = c("lightgrey", 'red'))        
  

combinedT_S2$Identity <- as.character(combinedT_S2$seurat_clusters)
combinedT_S2$Identity[!combinedT_S2$Identity %in% c('1')]<- 'Th2'
combinedT_S2$Identity[combinedT_S2$Identity %in% c('1')]<- 'CD4+CXCR5+T'
Idents(combinedT_S2) <- combinedT_S2$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(combinedT_S2, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)

combinedT_S3 <- combinedT[, combinedT$common_snn_res.0.6 == '1']
DefaultAssay(combinedT_S3) <- 'common'
combinedT_S3 <- FindNeighbors(combinedT_S3, reduction = "harmony", dims = 1:15)


combinedT_S3 <- FindClusters(combinedT_S3, resolution = 0.2, verbose = F)
options(repr.plot.width=6, repr.plot.height=6)
 DimPlot(combinedT_S3, label = T)

combinedT_S3$Identity <- as.character(combinedT_S3$seurat_clusters)
combinedT_S3$Identity[!combinedT_S3$Identity %in% c('2')]<- 'CD4+TN'
combinedT_S3$Identity[combinedT_S3$Identity %in% c('2')]<- 'CD4+CXCR5+T'

Idents(combinedT_S3) <- combinedT_S3$Identity
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(combinedT_S3, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)

combinedT$Identity <- as.character(combinedT$common_snn_res.0.6)
combinedT$Identity[combinedT$Identity %in% c('5')]<- 'Th1'
combinedT$Identity[combinedT$Identity %in% c('8')]<- 'Th2'
combinedT$Identity[combinedT$Identity %in% c('2')]<- 'CD4+CTL1'
combinedT$Identity[combinedT$Identity %in% c('11')]<- 'MKI67+T'

# combinedT$Identity[combinedT$Identity %in% c('5')]<- 'CD4+CXCR5+T'
combinedT$Identity[combinedT$Identity %in% c('0')]<- 'CD4+TCM'
combinedT$Identity[combinedT$Identity %in% c('12')]<- 'CD4+TN'
combinedT$Identity[combinedT$Identity %in% c('3')]<- 'Treg'
# combinedT$Identity[combinedT$Identity %in% c('10')]<- 'ILC2'
combinedT$Identity[combinedT$Identity %in% c('10')]<- 'CD4+CTL2'
combinedT$Identity[combinedT$Identity %in% c('4','9')]<- 'Th17'
Idents(combinedT) <- combinedT$Identity
combinedT <- Idents_proj(combinedT_S1, combinedT)
combinedT <- Idents_proj(combinedT_S2, combinedT)
combinedT <- Idents_proj(combinedT_S3, combinedT)
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(combinedT, reduction = 'umap',label = T, label.size = 20/.pt)+
scale_color_manual(values = color_used)



unique(combinedT$batch)



CD4T_metadata <- 
read.csv('/data1/02.private/dengyj/analysis/Hif2a/public_data/Hs_asthma_SI_2023/CD4T_metadata.csv',row.names = 1)



cluster_colors <- c('CD4+TN'='#F57F21','CD4+TCM'='#0FAA07','CD4+CXCR5+T'='#909ee7',#'#E41E25',
                    'Th2'='#89CC9C',
'ILC2'='#80B4E0','Treg'='#C35338','Th1'='#E05C6E','CD4+CTL1'='#22A3FF',
'CD4+CTL2'='#DCD44E','undefined'='#919191', 'CD4+TMC'='#CE999C','MKI67+T' = '#ed2126', 
                    'Th17' = '#c99c5d', 'SLPI+T'='#db72b4', 'S1PR1+Th2'='#70c6cd'
                   )





options(repr.plot.width=6, repr.plot.height=6)
DimPlot(combinedT, reduction = 'umap',label = T, label.size = 20/.pt, pt.size = 1, repel=T)+
theme(panel.border=element_rect(colour = "black", fill =NA),
          panel.background = element_rect(colour = NA, fill ='white'), 
      axis.title = element_text(size = 15),axis.text = element_text(size = 15),
      legend.title = element_text(size = 15),legend.text = element_text(size = 15),
      #legend.position = c(0.575, 0.9),
      #legend.direction = 'horizontal',
     plot.title = element_text(size = 15, face='italic', hjust=0.5))+
guides(color = F)+
scale_color_manual(values = c('CD4+TN'='#F57F21','CD4+TCM'='#0FAA07','CD4+CXCR5+T'='#909ee7',#'#E41E25',
                    'Th2'='#89CC9C',
                    'ILC2'='#80B4E0','Treg'='#C35338','Th1'='#E05C6E','CD4+CTL1'='#22A3FF',
                    'CD4+CTL2'='#DCD44E','undefined'='#919191', 'CD4+TMC'='#CE999C','MKI67+T' = '#ed2126', 
                    'Th17' = '#c99c5d', 'SLPI+T'='#db72b4', 'S1PR1+Th2'='#70c6cd'))

filtered_combinedT <- combinedT[, !Idents(combinedT) %in% c('ILC2')]



options(repr.plot.width=5, repr.plot.height=5)
p <- DimPlot(filtered_combinedT, reduction = 'umap',label = T, label.size = 17.5/.pt, pt.size = 1, repel=T)+
theme(panel.border=element_rect(colour = "black", fill =NA),
          panel.background = element_rect(colour = NA, fill ='white'), 
      axis.title = element_text(size = 17.5),axis.text = element_text(size = 17.5),
      legend.title = element_text(size = 17.5),legend.text = element_text(size = 17.5),
      #legend.position = c(0.575, 0.9),
      #legend.direction = 'horizontal',
     plot.title = element_text(size = 17.5, face='italic', hjust=0.5))+
guides(color = F)+
scale_color_manual(values = c('CD4+TN'='#F57F21','CD4+TCM'='#0FAA07','CD4+CXCR5+T'='#909ee7',#'#E41E25',
                    'Th2'='#89CC9C',
                    'ILC2'='#80B4E0','Treg'='#C35338','Th1'='#E05C6E','CD4+CTL1'='#22A3FF',
                    'CD4+CTL2'='#DCD44E','undefined'='#919191', 'CD4+TMC'='#CE999C','MKI67+T' = '#ed2126', 
                    'Th17' = '#c99c5d', 'SLPI+T'='#db72b4', 'S1PR1+Th2'='#70c6cd'))
p$layers[[2]]$data$ident <- gsub('CD4\\+', '', p$layers[[2]]$data$ident)
p$layers[[2]]$data$ident[p$layers[[2]]$data$ident=='TN'] <- 'Naive'
p


saveRDS(combinedT, file = 'combinedT.rds')
saveRDS(filtered_combinedT, file = 'filtered_combinedT.rds')