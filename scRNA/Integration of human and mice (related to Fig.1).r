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


