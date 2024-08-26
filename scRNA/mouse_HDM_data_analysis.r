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


####人鼠整合

library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)
library(openxlsx)
library(future)
library(org.Mm.eg.db)
library(SCopeLoomR)
library(plyr)

library(ComplexHeatmap)

library(patchwork)

library(ggforce)

source('/data1/02.private/dengyj/analysis/mycode/plot_helper.r')

this_study_Th2 <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/basic_analysis/Hif2a_Th2_2.rds')

#this_study_Th2<- this_study_Th2[rowSums(this_study_Th2$RNA@data)>0, ]

this_study_Th2$model <- 'OVA'
this_study_Th2$source <- 'this_study'



Mm_public_Th2 <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/public_data/Mm_combined/Th2.rds')



Hs_Th2 <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/public_data/Hs_combined/combined_Th2.rds')

Mm_common_features <- intersect(rownames(this_study_Th2), rownames(Mm_public_Th2$common))

#Mm_Th2 <- merge(this_study_Th2[, !this_study_Th2$seurat_clusters == '5'], Mm_public_Th2)
Mm_Th2 <- merge(this_study_Th2, Mm_public_Th2)
Mm_Th2[['common']] <- CreateAssayObject(counts = Mm_Th2$RNA@counts[Mm_common_features, ])
Mm_Th2$common@data <- Mm_Th2$RNA@data[rownames(Mm_Th2$common), ]



genes_tsv <- read.table('/data1/01.project/Hif2a//ABFC20200274-060709（ctm_hif2a）//RNA/ABFC20200274_10X-scRNA_result/2.Cellranger/W1-10X5/filtered_feature_bc_matrix/features.tsv.gz')
colnames(genes_tsv) <- c('ensembl', 'symbol')

genes_tsv$symbol <- make.unique(genes_tsv$symbol)

Hs_to_Mm <- 
data.table::fread('/data1/02.private/dengyj/analysis/database/specie_conversion/mart_export_ensemb_111_GRCh38.p14.txt')
Hs_to_Mm <- data.frame(Hs_to_Mm)

Hs_to_Mm <- Hs_to_Mm[!is.na(Hs_to_Mm$`Mouse.orthology.confidence..0.low..1.high.`), ]
#Hs_to_Mm <- Hs_to_Mm[Hs_to_Mm$`Mouse.orthology.confidence..0.low..1.high.` == 1, ]

exp <- Hs_Th2$common@data

exp_filtered <- exp[rownames(exp) %in% Hs_to_Mm$Gene.name, ]
#exp_filtered <- exp_filtered[rowSums(exp_filtered) > 0, ]

new_names <- lapply(rownames(exp_filtered), function(x){
    Mm_ensembl <- unique(Hs_to_Mm[Hs_to_Mm$Gene.name %in% x, 'Mouse.gene.stable.ID'])
    if(length(Mm_ensembl) > 0){
        ###注意，会出现一个人源基因匹配多个小鼠基因的情况
        #Mm_symbol <- unique(genes_tsv[genes_tsv$ensembl %in% Mm_ensembl, 'symbol'])
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

all(rownames(Mm_Th2$common@data) %in% genes_tsv$symbol)

Mm_Th2_ensembl_exp <- Mm_Th2$common@data

Mm_Th2_ensembl_id <- 
get_id_changed(rownames(Mm_Th2$common@data), mapping = genes_tsv, orig_var = 'symbol', new_var = 'ensembl')

rownames(Mm_Th2_ensembl_exp) <- Mm_Th2_ensembl_id

Ms_symbol_ensembl <- data.frame(Ms_symbol = rownames(Mm_Th2$common@data), Mm_ensembl = 
                               Mm_Th2_ensembl_id)

sum(duplicated(rownames(Mm_Th2_ensembl_exp)))

Mm_Th2_ensembl_exp[1:10,1:10]



HsMm_common_features <- intersect(rownames(exp_filtered_Mm), rownames(Mm_Th2_ensembl_exp))

length(HsMm_common_features)

Hs_mapping_Mm <- Hs_mapping_Mm[Hs_mapping_Mm$Mm_ensembl %in% HsMm_common_features, ]
Ms_symbol_ensembl <- Ms_symbol_ensembl[Ms_symbol_ensembl$Mm_ensembl %in% HsMm_common_features, ]

mapping_df <- merge(Hs_mapping_Mm, Ms_symbol_ensembl)

# tmp <- fread('/data1/02.private/dengyj/analysis/database/GSEA/Mouse_Gene_Symbol_Remapping_MSigDB.v7.0.chip')

# head(tmp)

new_Hs_Th2 <- CreateSeuratObject(counts = exp_filtered_Mm[HsMm_common_features, ], meta.data = Hs_Th2@meta.data)
new_Hs_Th2$species <- 'human'
new_Hs_Th2$batch_to_reduce <- new_Hs_Th2$source
new_Hs_Th2$batch_to_reduce[new_Hs_Th2$batch_to_reduce %in% c('SI', 'NatMed_asthma')] <- 'Asthma'
new_Hs_Th2$orig_idents <- as.character(Idents(Hs_Th2))

Mm_Th2_ensembl <- CreateSeuratObject(counts = Mm_Th2_ensembl_exp[HsMm_common_features, ], meta.data = Mm_Th2@meta.data)
Mm_Th2_ensembl$species <- 'mouse'
Mm_Th2_ensembl$batch_to_reduce <- Mm_Th2_ensembl$source

logi <- Mm_Th2_ensembl$batch_to_reduce == 'this_study'
Mm_Th2_ensembl$batch_to_reduce[logi] <- Mm_Th2_ensembl$batch[logi]

Mm_Th2_ensembl$orig_idents <- as.character(Mm_Th2_ensembl$seurat_clusters)

HsMm_Th2 <- merge(Mm_Th2_ensembl, new_Hs_Th2)

dim(new_Hs_Th2)

HsMm_Th2_symbol_mat <- HsMm_Th2$RNA@data

rownames(HsMm_Th2_symbol_mat) <- mapping_df[match(rownames(HsMm_Th2_symbol_mat), mapping_df$Mm_ensembl), 'Ms_symbol']

HsMm_Th2[['symbol']] <- CreateAssayObject(counts = HsMm_Th2_symbol_mat)

DefaultAssay(HsMm_Th2) <- 'RNA'
HsMm_Th2_list <- SplitObject(HsMm_Th2, split.by = 'batch_to_reduce')
HsMm_Th2_list <- lapply(HsMm_Th2_list, function(obj){
    obj <- FindVariableFeatures(obj)
    obj
})

dim(HsMm_Th2[['symbol']])

anchors_features <- SelectIntegrationFeatures(HsMm_Th2_list, nfeatures = 2500)#OVA_Th2@assays$RNA@var.features

all(anchors_features %in% rownames(HsMm_Th2$RNA))

table(HsMm_Th2$batch_to_reduce)

options(future.globals.maxSize = 10 * 1024^3)
anchors <- FindIntegrationAnchors(object.list = HsMm_Th2_list, dims = 1:20, verbose = F,
                                  anchor.features = anchors_features, k.filter = 150)
HsMm_Th2 <- IntegrateData(anchorset = anchors, dims = 1:20, verbose = F)

DefaultAssay(HsMm_Th2) <- 'integrated'
HsMm_Th2 <- ScaleData(HsMm_Th2, verbose  = F)#, vars.to.regress = names(Go_list))
HsMm_Th2 <- RunPCA(HsMm_Th2, npcs = 50, verbose = T)

options(repr.plot.width=6, repr.plot.height=6)
ElbowPlot(HsMm_Th2, ndims = 50)

HsMm_Th2 <- RunUMAP(HsMm_Th2, reduction = "pca", dims = setdiff(1:10, c()), verbose = F)

options(repr.plot.width=8, repr.plot.height=8)
HsMm_Th2 <- AddMetaData(HsMm_Th2, this_study_Th2$seurat_clusters, col.name = 'orig_id')
DimPlot(HsMm_Th2[, HsMm_Th2$source %in% c('this_study')], reduction = "umap", 
        label = T, group.by = 'orig_id')



options(repr.plot.width=8, repr.plot.height=6)
DimPlot(HsMm_Th2, reduction = "umap", label = T, group.by = 'batch')

options(repr.plot.width=8, repr.plot.height=6)
DimPlot(HsMm_Th2, reduction = "umap", label = T, group.by = 'species')

DefaultAssay(HsMm_Th2) <- 'symbol'
options(repr.plot.width=16, repr.plot.height=16)

FeaturePlot(HsMm_Th2, reduction = "umap", features = c('Cxcr6', 'Cxcr5','Bcl6', 'Pdcd1'))

DefaultAssay(HsMm_Th2) <- 'symbol'
options(repr.plot.width=16, repr.plot.height=16)

FeaturePlot(HsMm_Th2, reduction = "umap", features = c('Il1rl1', 'Epas1', 'Il13', 'Il5'))

DefaultAssay(HsMm_Th2) <- 'symbol'
options(repr.plot.width=16, repr.plot.height=16)

FeaturePlot(HsMm_Th2, reduction = "umap", features = c('Il17rb', 'Plin2', 'Pparg', 'Cxcr6'))

DefaultAssay(HsMm_Th2) <- 'symbol'
options(repr.plot.width=16, repr.plot.height=16)

FeaturePlot(HsMm_Th2, reduction = "umap", features = c('S100a4', 'Il18r1', 'Ctsw', 'Itage'))

DefaultAssay(HsMm_Th2) <- 'symbol'
options(repr.plot.width=16, repr.plot.height=16)

FeaturePlot(HsMm_Th2, reduction = "umap", features = c('Ctla4', 'Lag3', 'Foxp3', 'Cd27'))

options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(HsMm_Th2) <- 'symbol'
FeaturePlot(HsMm_Th2, features = c('Nr4a1'))

options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(HsMm_Th2) <- 'symbol'
FeaturePlot(HsMm_Th2, features = c('Isg15'))

DefaultAssay(HsMm_Th2) <- 'symbol'
options(repr.plot.width=16, repr.plot.height=16)

FeaturePlot(HsMm_Th2, reduction = "umap", features = c('Tcf7', 'Tox', 'Tox2', 'Btla'))

options(repr.plot.width=8, repr.plot.height=8)
FeaturePlot(HsMm_Th2, reduction = "umap", features = c('Slamf6'))

options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(HsMm_Th2) <- 'symbol'
FeaturePlot(HsMm_Th2, features = c('Tcf7'))

options(repr.plot.width=8, repr.plot.height=8)
FeaturePlot(HsMm_Th2, reduction = "umap", features = c('S1pr1'))

colors <- c('#F57F21', '#919191','#0FAA07', '#C35338', '#89CC9C', '#80B4E0',
 '#DCD44E', '#E05C6E', '#64446C', '#22A3FF', '#E41E25','#CE999C',  
 '#136433', '#9C56A3', '#FE3DF6', '#9F8CFF')



cluster_colors <- 
c("Th2"="#E5D2DD","T_memory"="#53A85F","Th17"="#F1BB72","Treg"="#F3B1A0",'ST2+Treg'='#859472',
  "Alveolar_Mac"="#D6E7A3",
"IFN_T"="#57C3F3","Th1"="#476D87","Cd8_T"="#E95C59","Trans_mac"="#E59CC4","T_naive"="#AB3282",
"Neutrophil"="#23452F","Bcell"="#BD956A","Inflammatory_Mac"="#8C549C","Mki67_T"="#585658",
"Eosinophil"="#9FA3A8","Interstitial_mac"="#E0D4CA","Epi"="#5F3D69","cDC"="#C5DEBA",
  "γδT"="#58A4C3",
"Mast"="#E4C755","pDC"="#F7F398","17"="#AA9A59","8"="#E63863","other" = "gray80", 
 'Tcell'='#4284c9', 
# "0"="#C61584","1"="#AAD0BA","2"="#5CA553","3"="#E78C5D","4"="#66AAE2","5"="#F9EE40"
 "0"="#E3D6E0",
  "1"="#5B8FAE","2"="#C667AA","3"="#E8646D","4"="#65AB66","5"="#B79EE8", 
                   'T_naive'= '#dcc442')

DefaultAssay(HsMm_Th2) <- 'integrated'
HsMm_Th2 <- FindNeighbors(HsMm_Th2, dims = 1:10)

HsMm_Th2 <- FindClusters(HsMm_Th2, resolution = seq(0.2, 1, 0.05), verbose = F)

options(repr.plot.width=8, repr.plot.height=8)
for(i in seq(0.2, 1, 0.05)){
    res <- paste0('integrated_snn_res.', i)
    p <- DimPlot(HsMm_Th2, reduction = 'umap',label = T, label.size = 20/.pt, 
                 group.by = res)+
    labs(title = res)

    print(p)

}



options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(HsMm_Th2) <- 'symbol'
FeaturePlot(HsMm_Th2[, HsMm_Th2$species =='human'], features = c('Il1rl1', 'Il13', 'Gata3','Il5'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(HsMm_Th2) <- 'symbol'
FeaturePlot(HsMm_Th2[, HsMm_Th2$species =='mouse'], features = c('Il1rl1', 'Il13', 'Gata3','Il5'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(HsMm_Th2) <- 'symbol'
FeaturePlot(HsMm_Th2, features = c('Il10', 'Ctla4', 'Cd27', 'Foxp1'))

options(repr.plot.width=8, repr.plot.height=8)
DefaultAssay(HsMm_Th2) <- 'symbol'
FeaturePlot(HsMm_Th2, features = c('Il2ra'))

options(repr.plot.width=16, repr.plot.height=16)
DefaultAssay(HsMm_Th2) <- 'symbol'
FeaturePlot(HsMm_Th2, features = c('Isg15', 'Mx1', 'Usp18', 'Dgat1'))



table(filtered_HsMm_Th2$species)

Idents(HsMm_Th2) <- HsMm_Th2$integrated_snn_res.0.4
options(repr.plot.width=8, repr.plot.height=8)
DimPlot(HsMm_Th2, reduction = 'umap',label = T, label.size = 20/.pt)

options(repr.plot.width=8, repr.plot.height=8)
HsMm_Th2$seu_ID <- as.character(HsMm_Th2$integrated_snn_res.0.4)
HsMm_Th2$seu_ID[HsMm_Th2$integrated_snn_res.0.4 %in% c('0','2')] <- '3'
HsMm_Th2$seu_ID[HsMm_Th2$integrated_snn_res.0.4 %in% c('1')] <- '1'
HsMm_Th2$seu_ID[HsMm_Th2$integrated_snn_res.0.4 %in% c('5','6','7')] <- '0'

HsMm_Th2$seu_ID[HsMm_Th2$integrated_snn_res.0.4 %in% c('3')] <- '2'
HsMm_Th2$seu_ID[HsMm_Th2$integrated_snn_res.0.4 %in% c('4')] <- '4'
HsMm_Th2$seu_ID[HsMm_Th2$integrated_snn_res.0.4 %in% c('8')] <- '5'
DimPlot(HsMm_Th2, reduction = 'umap',label = T, label.size = 20/.pt, 
                 group.by = 'seu_ID')+
scale_color_manual(values = cluster_colors)

HsMm_Th2$seu_ID <- factor(HsMm_Th2$seu_ID, levels = c('0','1','2','3','4','5'))



table(filtered_HsMm_Th2$species)

t(filtered_HsMm_Th2$source)





HsMm_Th2$source2 <- ''
HsMm_Th2$source2[HsMm_Th2$source =='this_study'] <- 'this_study'
HsMm_Th2$source2[HsMm_Th2$source !='this_study' & HsMm_Th2$species == 'mouse'] <- 'public (mouse)'
HsMm_Th2$source2[HsMm_Th2$source !='this_study' & HsMm_Th2$species == 'human'] <- 'public (human)'
HsMm_Th2$source2 <- factor(HsMm_Th2$source2, levels = c('this_study', 'public (mouse)', 
                                                       'public (human)'))

filtered_HsMm_Th2 <- HsMm_Th2[, !HsMm_Th2$source2 %in% c('this_study')]

cluster_colors <- 
c("Th2"="#E5D2DD","T_memory"="#53A85F","Th17"="#F1BB72","Treg"="#F3B1A0",'ST2+Treg'='#859472',
  "Alveolar_Mac"="#D6E7A3",
"IFN_T"="#57C3F3","Th1"="#476D87","Cd8_T"="#E95C59","Trans_mac"="#E59CC4","T_naive"="#AB3282",
"Neutrophil"="#23452F","Bcell"="#BD956A","Inflammatory_Mac"="#8C549C","Mki67_T"="#585658",
"Eosinophil"="#9FA3A8","Interstitial_mac"="#E0D4CA","Epi"="#5F3D69","cDC"="#C5DEBA",
  "γδT"="#58A4C3",
"Mast"="#E4C755","pDC"="#F7F398","17"="#AA9A59","8"="#E63863","other" = "gray80", 
 'Tcell'='#4284c9', 

  '0'='#91D0D4', '1'='#587DB2', '2'='#42827E',
                   
 '3'='#A35592', '4'='#E9A342','5'='#DBACC2',
                   'T_naive'= '#dcc442')

plot_df <- cbind(Embeddings(filtered_HsMm_Th2, reduction = 'umap'), data.frame(identity=Idents(filtered_HsMm_Th2)))
options(repr.plot.width=5, repr.plot.height=5)
umap_p_idents <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = identity))+
geom_point(size=1)+
  theme(legend.position = c(0.17,0.1),
        panel.background=element_rect(colour = NA, fill = 'white'),
        panel.border=element_rect(colour = 'black', fill = NA),
        axis.title=element_text(size = 20),axis.ticks=element_blank(),
        legend.title=element_blank(),
        legend.background=element_blank(),
        axis.text=element_blank(),axis.line=element_blank(), 
        legend.text = element_text(size = 20))+#ylim(-6,6)+xlim(-6,7.5)+
  guides(colour = guide_legend(override.aes = list(size=7), ncol = 1)) +
  scale_color_manual(values=cluster_colors, label = 
                    function(x){
                        paste0('C', x)
                    })+
guides(color = guide_legend(nrow = 3, override.aes = list(size = 5)))
umap_p_idents


###########Correlation between EPAS1 and PATH2 signature


library(Seurat)
library(ggplot2)
library(VISION)

source('/data1/02.private/dengyj/analysis/Hif2a/plot/support/plot_helper.r')

Hs_public_Th2 <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/public_data/Hs_combined/combined_Th2.rds')
Hs_public_CD4T_es <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/public_data/Hs_combined/combinedT_es.rds')

Mm_public_Th2 <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/public_data/Mm_combined/Th2.rds')
Mm_public_CD4T_es <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/public_data/Mm_combined/combined_CD4T_es.rds')






Hs_cor_reuslt_df <- lapply(seq(5, 50, 5), function(i){
    set.seed(1234)
    Hs_public_Th2_pools <- applyMicroClustering(Hs_public_Th2$RNA@data, cellsPerPartition = i, 
                                                latentSpace = Embeddings(Hs_public_Th2, 'harmony')[, 1:5])

    Hs_Th2_microcluster_id <- rep(names(Hs_public_Th2_pools), lengths(Hs_public_Th2_pools))
    names(Hs_Th2_microcluster_id) <- unlist(Hs_public_Th2_pools)
    Hs_public_Th2 <- AddMetaData(Hs_public_Th2, Hs_Th2_microcluster_id, col.name = 'Hs_Th2_microcluster_id')
    Hs_public_Th2 <- AddMetaData(Hs_public_Th2, t(Hs_public_CD4T_es), col.name = 'ST2_hi_Th2')

    Hs_data <- t(cbind(gene = expm1(Hs_public_Th2$RNA@data['EPAS1',]), 
                     FetchData(Hs_public_Th2, vars = c('ST2_hi_Th2'))))
    Hs_data_mean <- my_group_mean(Hs_data, factor(Hs_public_Th2$Hs_Th2_microcluster_id))
    Hs_cor_result <- cor(t(Hs_data_mean), method = 'pearson')[1,2]
    df <- data.frame(iter = i, cor = Hs_cor_result)
    df
    #print(paste0('iter: ', i,' cor: ', Hs_cor_result))   
})


Hs_cor_reuslt_df <- do.call(rbind, Hs_cor_reuslt_df)


Mm_cor_reuslt_df <- lapply(seq(5, 50, 5), function(i){
    set.seed(1234)
    Mm_public_Th2_pools <- applyMicroClustering(Mm_public_Th2$RNA@data, cellsPerPartition = i, 
                                                latentSpace = Embeddings(Mm_public_Th2, 'harmony')[, 1:10])

    Mm_Th2_microcluster_id <- rep(names(Mm_public_Th2_pools), lengths(Mm_public_Th2_pools))
    names(Mm_Th2_microcluster_id) <- unlist(Mm_public_Th2_pools)
    Mm_public_Th2 <- AddMetaData(Mm_public_Th2, Mm_Th2_microcluster_id, col.name = 'Mm_Th2_microcluster_id')
    Mm_public_Th2 <- AddMetaData(Mm_public_Th2, t(Mm_public_CD4T_es), col.name = 'ST2_hi_Th2')

    Mm_data <- t(cbind(gene = expm1(Mm_public_Th2$RNA@data['Epas1',]), 
                     FetchData(Mm_public_Th2, vars = c('ST2_hi_Th2'))))
    Mm_data_mean <- my_group_mean(Mm_data, factor(Mm_public_Th2$Mm_Th2_microcluster_id))
    Mm_cor_result <- cor(t(Mm_data_mean), method = 'pearson')[1,2]
    df <- data.frame(iter = i, cor = Mm_cor_result)
    df
    #print(paste0('iter: ', i,' cor: ', Mm_cor_result))   
})


Mm_cor_reuslt_df <- do.call(rbind, Mm_cor_reuslt_df)

Hs_cor_reuslt_df

Mm_cor_reuslt_df





set.seed(1234)

Hs_public_Th2_pools <- applyMicroClustering(Hs_public_Th2$RNA@data, cellsPerPartition = 45, 
                                            latentSpace = Embeddings(Hs_public_Th2, 'harmony')[, 1:5])

Hs_Th2_microcluster_id <- rep(names(Hs_public_Th2_pools), lengths(Hs_public_Th2_pools))
names(Hs_Th2_microcluster_id) <- unlist(Hs_public_Th2_pools)
Hs_public_Th2 <- AddMetaData(Hs_public_Th2, Hs_Th2_microcluster_id, col.name = 'Hs_Th2_microcluster_id')
Hs_public_Th2 <- AddMetaData(Hs_public_Th2, t(Hs_public_CD4T_es), col.name = 'ST2_hi_Th2')

Hs_data <- t(cbind(gene = expm1(Hs_public_Th2$RNA@data['EPAS1',]), 
                 FetchData(Hs_public_Th2, vars = c('ST2_hi_Th2'))))
Hs_data_mean <- my_group_mean(Hs_data, factor(Hs_public_Th2$Hs_Th2_microcluster_id))

Hs_plot_df <- data.frame(t(Hs_data_mean))



set.seed(1234)
Mm_public_Th2_pools <- applyMicroClustering(Mm_public_Th2$RNA@data, cellsPerPartition = 45, 
                                            latentSpace = Embeddings(Mm_public_Th2, 'harmony')[, 1:10])

Mm_Th2_microcluster_id <- rep(names(Mm_public_Th2_pools), lengths(Mm_public_Th2_pools))
names(Mm_Th2_microcluster_id) <- unlist(Mm_public_Th2_pools)
Mm_public_Th2 <- AddMetaData(Mm_public_Th2, Mm_Th2_microcluster_id, col.name = 'Mm_Th2_microcluster_id')
Mm_public_Th2 <- AddMetaData(Mm_public_Th2, t(Mm_public_CD4T_es), col.name = 'ST2_hi_Th2')

Mm_data <- t(cbind(gene = expm1(Mm_public_Th2$RNA@data['Epas1',]), 
                 FetchData(Mm_public_Th2, vars = c('ST2_hi_Th2'))))
Mm_data_mean <- my_group_mean(Mm_data, factor(Mm_public_Th2$Mm_Th2_microcluster_id))
Mm_plot_df <- data.frame(t(Mm_data_mean))

options(repr.plot.width = 5, repr.plot.height = 4.75)
ggplot(Hs_plot_df, aes(gene, ST2_hi_Th2))+
geom_point(color ='#336caf', size = 5, shape = 19)+
geom_smooth(method = 'lm', se = T, color ='red')+
ggpubr::stat_cor(method = "pearson",size=20/.pt, label.y.npc = 0.2, label.x.npc = 0.22, color = 'red')+
labs(y='PATH2 score', x = 'Expression level of EPAS1')+
    theme(axis.ticks = element_blank(), 
          plot.margin = margin(r = 0.45, t = 0.1,unit = 'cm'),
          axis.line = element_blank(),
         axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, face = 'italic'),
          plot.title = element_text(size = 20,hjust = 0.5),
          axis.text = element_text(size = 20),
          legend.justification = 0.5,
        legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          #legend.position = c(0.9, 0.2),
          panel.border = element_rect(color ='black', fill =NA),
          panel.background= element_rect(fill ='white', color =NA),
          
         )+
guides(color = guide_legend(override.aes = list(size=5)))+
labs(title = 'Pearson correlation (human)')
ggsave('Epas1_Path2_cor_human.svg', width = 5, height = 4.75)

?ggpubr::stat_cor

getwd()

options(repr.plot.width = 5, repr.plot.height = 4.75)
ggplot(Mm_plot_df, aes(gene, ST2_hi_Th2))+
geom_point(color ='#4baf33', size = 5, shape = 17)+
geom_smooth(method = 'lm', se = T, color ='red')+
ggpubr::stat_cor(method = "pearson",size=20/.pt, label.y.npc = 0.06, label.x.npc = 0.25, color = 'red')+
labs(y='PATH2 score', x = 'Expression level of Epas1')+
    theme(axis.ticks = element_blank(),
          axis.line = element_blank(),
         axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, face = 'italic'),
          plot.title = element_text(size = 20,hjust = 0.5),
          axis.text = element_text(size = 20),
          legend.justification = 0.5,
        legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          #legend.position = c(0.9, 0.2),
          panel.border = element_rect(color ='black', fill =NA),
          panel.background= element_rect(fill ='white', color =NA),
          
         )+
guides(color = guide_legend(override.aes = list(size=5)))+
labs(title = 'Pearson correlation (mouse)')
ggsave('Epas1_Path2_cor_mouse.svg', width = 5, height = 4.75)

getwd()

40*4.75/5

getwd()

getwd()

save(Hs_plot_df, Mm_plot_df, file = 'Epas1_Path2_cor.RData')

save(Hs_public_Th2_pools, Mm_public_Th2_pools, file = 'public_Th2_pools.RData')

