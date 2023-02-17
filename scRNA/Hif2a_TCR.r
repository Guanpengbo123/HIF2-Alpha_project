#####TCR diversity between Th2 subsets and circle packing plot
library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(treemapify)
library(dplyr)
library(plyr)
library(ggpubr)
source('/data1/02.private/dengyj/analysis/Hif2a/TCR/TCR_caculation.R')
library(patchwork)
library(rlang)
library(reshape2)
library(openxlsx)
library(packcircles)
library(ggforce)


source('/data1/02.private/dengyj/analysis/mycode/GeomLabelCircle.R')

setdiff_clu <- c('γδT', 'Eosinophil', 'Interstitial_mac','Alveolar_Mac','Bcell', 'Inflammatory_Mac', 
 'Neutrophil','Epi', 'Trans_mac','pDC','cDC','Mast', 'Th2')

table(hif_step_2$batch)

Hif2a_Th2_2 <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/basic_analysis/Hif2a_Th2_2.rds')

Idents(Hif2a_Th2_2) <- Hif2a_Th2_2$identity

meta <- Hif2a_Th2_2[[]]

tcr_dir1 <- '/data3/Hif2a/ABFC20200274-161728/ABFC20200274_10X-TCR_result/2.Cellranger/'
tcr_dir1_files <- dir(tcr_dir1)
tcr_dir1_files <- tcr_dir1_files[grep('TCR', tcr_dir1_files)]

tcr_dir1 <- '/data3/Hif2a/ABFC20200274-161728/ABFC20200274_10X-TCR_result/2.Cellranger/'
tcr_dir1_files <- dir(tcr_dir1)
tcr_dir1_files <- tcr_dir1_files[grep('TCR', tcr_dir1_files)]
#tcr_dir1_files <- tcr_dir1_files[!grepl('I1', tcr_dir1_files)]
tcr_raw_data1 <- lapply(tcr_dir1_files, function(x){
    tmp_dir <- paste0(tcr_dir1, x, '/filtered_contig_annotations.csv')
    df <- read.csv(tmp_dir)
    sample <- gsub('TCR', '5', x)
    df$barcode <- paste0(df$barcode, '_', sample)
    df <- df[df$barcode %in% rownames(meta), ]###所有细胞都在meta
    df <- df[grepl('true', df$high_confidence, ignore.case = T), ]
    df <- df[grepl('true', df$productive, ignore.case = T), ]
    if(nrow(df) > 0){
        df$Sample <- meta[match(df$barcode, rownames(meta)), 'batch']
        df$identity <- meta[match(df$barcode, rownames(meta)), 'identity']
        #df <- df[!df$identity %in% setdiff_clu, ]
        df$Source <- meta[match(df$barcode, rownames(meta)), 'Source']
        df$raw_clonotype_id <- paste0(df$Sample, '_', df$raw_clonotype_id)
        return(df)        
    }

})


tcr_dir2 <- '/data3/Hif2a/W1_W2_K2_K3_TCR/TCR/3.Cellranger/'
tcr_dir2_files <- dir(tcr_dir2)
tcr_dir2_files <- tcr_dir2_files[grep('TCR', tcr_dir2_files)]
#tcr_dir2_files <- tcr_dir2_files[!grepl('I1', tcr_dir2_files)]
tcr_raw_data2 <- lapply(tcr_dir2_files, function(x){
    tmp_dir <- paste0(tcr_dir2, x, '/', x, '/outs/filtered_contig_annotations.csv')
    df <- read.csv(tmp_dir)
    sample <- gsub('TCR', '5', x)
    df$barcode <- paste0(df$barcode, '_', sample)
    df <- df[df$barcode %in% rownames(meta), ]
    df <- df[grepl('true', df$high_confidence, ignore.case = T), ]
    df <- df[grepl('true', df$productive, ignore.case = T), ]    
    if(nrow(df) > 0){
        df$Sample <- meta[match(df$barcode, rownames(meta)), 'batch']
        df$identity <- meta[match(df$barcode, rownames(meta)), 'identity']
        #df <- df[!df$identity %in% setdiff_clu, ]
        df$Source <- meta[match(df$barcode, rownames(meta)), 'Source']
        df$raw_clonotype_id <- paste0(df$Sample, '_', df$raw_clonotype_id)
        return(df) 
    }
    
})

tcr_raw <- list(data = c(tcr_raw_data1, tcr_raw_data2), meta = meta)
#tcr_raw <- list(data = c(tcr_raw_data1, tcr_raw_data2), meta = meta)
names(tcr_raw$data) <- c(tcr_dir1_files, tcr_dir2_files)
tcr_raw$data <- tcr_raw$data[!sapply(tcr_raw$data, is.null)]


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




div_by_cells_idents <- lapply(tcr_raw$data, 
                           function(x) get_div_by_cells(df = x, group_var = 'identity', method = 'div')
                          )
div_by_cells_idents <- do.call(rbind, div_by_cells_idents)
div_by_cells_idents <- div_by_cells_idents[div_by_cells_idents$identity %in% 
                                           c('0', '1', '2', '3', '4', '5'), ]
div_by_cells_idents$source <- ifelse(test = grepl('^K', rownames(div_by_cells_idents)), 
                                    yes = 'KO', no = 'WT')

###############
options(repr.plot.width=8, repr.plot.height=8)

ggplot(div_by_cells_idents,#[div_by_cells_idents$identity %in% clusters_used, ], 
       aes(identity, Diversity, fill = identity))+ geom_boxplot(outlier.alpha = 0)+
stat_compare_means( method="kruskal.test", hjust =-0.5, size = 15/.pt, vjust = -0.5)
geom_point(mapping = aes(color =source), position = position_dodge(0.9))+
scale_fill_manual(values = cluster_colors)+
labs(y = 'TCR diversity per cella', fill = 'identity', x = NULL)





############circle packing
tcr_df <- do.call(rbind, tcr_raw$data)
tcr_df <- tcr_df[!duplicated(tcr_df$barcode),]

tcr_clone_clu_df <- data.frame(table(tcr_df$raw_clonotype_id, tcr_df$identity))
colnames(tcr_clone_clu_df) <- c('clonal_id', 'identity', 'count')
tcr_clone_clu_df <- tcr_clone_clu_df[tcr_clone_clu_df$count > 0, ]

Hif2a_Th2_2_umap_df <- cbind(FetchData(Hif2a_Th2_2, vars = c('identity')),
                             Embeddings(Hif2a_Th2_2, reduction = 'umap'))
Hif2a_Th2_2_umap_df <- ddply(Hif2a_Th2_2_umap_df, c('identity'),
                             summarise,UMAP_1=median(UMAP_1), UMAP_2=median(UMAP_2))

scale_FC <- 15
packing_df <- lapply(levels(tcr_clone_clu_df$identity), function(id){
    
    tmp_df <- tcr_clone_clu_df[tcr_clone_clu_df$identity %in% c(id), ]
    areas <- sort(tmp_df$count, decreasing = T)
    packing <- circleProgressiveLayout(areas, sizetype = 'area')
    
    packing$size <- packing$radius^2*pi
    median_loc <- c(median(packing$x), median(packing$y))
    new_median_loc <- Hif2a_Th2_2_umap_df[Hif2a_Th2_2_umap_df$identity %in% id, c('UMAP_1', 'UMAP_2')]
    new_packing <- lapply(1:nrow(packing), function(row){
        new_x<- unname((packing[row, 'x'] - median_loc[1])/scale_FC + new_median_loc[1])
        new_y<- unname((packing[row, 'y'] - median_loc[2])/scale_FC + new_median_loc[2])
        #new_size <- packing[row, 'size']/scale_FC^2
        new_df <- data.frame(x = new_x, y=new_y, 
                             size = packing[row, 'size'], #new_size= new_size, 
                            identity=id)
        new_df
    })
    new_packing <- do.call(rbind, new_packing)
    
    
})
names(packing_df) <- levels(tcr_clone_clu_df$identity)

new_packing_df <- do.call(rbind, packing_df)

options(repr.plot.width=8, repr.plot.height=8)
tmp_range <- range(new_packing_df$size)*0.3
tmp_range[1] <- tmp_range[1]*3.3
ggplot(data = new_packing_df, aes(x, y, size = size, fill = identity)) + 
 geom_point(shape=21, stroke=0.1, color='black') +
scale_size_continuous(range = tmp_range)+
  coord_equal() +
scale_fill_manual(values = cluster_colors)+
theme(panel.border=element_blank(),
          panel.background = element_blank(), 
      axis.ticks = element_blank(),
      axis.title = element_blank(),axis.text = element_blank(),
      legend.title = element_text(size = 15),legend.text = element_text(size = 15),

)+
guides(fill=guide_legend(override.aes = list(size=7)), 
       size= guide_legend(override.aes = list(fill='black')))+
labs(x=NULL,y=NULL, size = 'Clone size', fill ='Identity')


options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = packing_df[['1']], aes(x, y, size = size, fill = identity)) + 
 geom_point(shape=21, stroke=0.1, color='black') +
scale_size_area(max_size = max(packing_df[['1']]$size*0.95))+
  coord_equal() +
scale_fill_manual(values = cluster_colors)+
theme(panel.border=element_blank(),
          panel.background = element_blank(), 
      axis.ticks = element_blank(),
      axis.title = element_blank(),axis.text = element_blank(),
      legend.title = element_text(size = 15),legend.text = element_text(size = 15),

)+
guides(fill=guide_legend(override.aes = list(size=7)), 
       size= guide_legend(override.aes = list(fill='black')))+
labs(x=NULL,y=NULL, size = 'Clone size', fill ='Identity')
