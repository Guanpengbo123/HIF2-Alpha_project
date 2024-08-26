library(Seurat)
library(dplyr)

library(ggplot2)

library(openxlsx)
library(treemapify)

library(dplyr)

library(plyr)

library(ggpubr)
# library(cowplot)
source('/data1/02.private/dengyj/analysis/Hif2a/TCR/TCR_caculation.R')
library(patchwork)
library(rlang)
# library(ggalluvial)
library(reshape2)

library(openxlsx)

library(packcircles)

library(ggforce)

get_sem <- function(x){
    sem <- sd(x)/sqrt(length(x))
    val <- c((mean(x) -sem), (mean(x)+sem))
    names(val) <- c('ymin', 'ymax')
    return(val)
}

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

source('/data1/02.private/dengyj/analysis/mycode/GeomLabelCircle.R')

setdiff_clu <- c('γδT', 'Eosinophil', 'Interstitial_mac','Alveolar_Mac','Bcell', 'Inflammatory_Mac', 
 'Neutrophil','Epi', 'Trans_mac','pDC','cDC','Mast', 'Th2')
####Th2重新命名为了0，1，2，3，5，还命名为Th2的是亚聚类过程中的低质量细胞

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


tree_color <- c("#9E6BAB", "#CFE3C6", "#E3CEE4", "#F3746C", "#86382A", "#ABA300", "#D6D4EB"
, "#B9A96B", "#FF68A1", "#EAA944", "#7AAF93", "#7A8DBB", "#7673AE", "#396F68"
, "#ECAFCF", "#EACF68", "#F7DDD4", "#EBDBE4", "#66CEF6", "#F8F4A8", "#C35338"
, "#EF5276", "#A0D7C9", "#63B472", "#F9DBE5", "#0CB702", "#F48930", "#6B6A6B"
, "#27BDCF", "#F8BFAF", "#F5C785", "#DEEAB1", "#217CB8", "#31FEB3", "#74517B"
, "#588198", "#CAA57D", "#9C99C4", "#2D563D", "#FF77AB", "#9F8CFF", "#D5E7F7"
, "#22A3FF", "#00E8F7", "#BB4A94", "#69B4CE", "#C9BDB2")

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
                           function(x) get_div_by_cells(df = x, group_var = 'identity', method = 'pielou')
                          )
div_by_cells_idents <- do.call(rbind, div_by_cells_idents)
div_by_cells_idents <- div_by_cells_idents[div_by_cells_idents$identity %in% 
                                           c('0', '1', '2', '3', '4', '5'), ]
div_by_cells_idents$source <- ifelse(test = grepl('^K', rownames(div_by_cells_idents)), 
                                    yes = 'KO', no = 'WT')


options(repr.plot.width=6, repr.plot.height=6)

ggplot(div_by_cells_idents,#[div_by_cells_idents$identity %in% clusters_used, ], 
       aes(identity, Diversity, fill = identity))+ geom_boxplot(outlier.alpha = 0)+
stat_compare_means( method="kruskal.test", label.x.npc = 0.3, size = 15/.pt, vjust = -0.5)+
theme(axis.text.x = element_text(angle = 0, hjust=0.5, size = 14.5), 
     axis.text.y = element_text(size = 15), 
     axis.title.y = element_text(size = 15),
     legend.title = element_text(size = 20),legend.text = element_text(size = 15), 
      panel.background = element_rect(fill = 'white', color = NA), 
      panel.border = element_rect(color = 'black', fill = NA),
      plot.title = element_text(size = 15, hjust = 0.5)
     )+
scale_x_discrete(label = function(x){paste0('C', x)})+
geom_point(mapping = aes(color =source), position = position_dodge(0.9))+
guides(fill = guide_legend(override.aes = list(size=0.1), keywidth = 2, keyheight = 2), 
      color =guide_legend(override.aes = list(size=7)))+
scale_fill_manual(values = cluster_colors, label = function(x){
    x <- paste0('C', x)
})+
scale_color_manual(values = c('WT' = '#96d8dc', 'KO' = '#cabae4'), limits = c('WT', 'KO'))+
labs(y = "Pielou's evenness index", fill = 'identity', x = NULL)
ggsave('Th2_subsets_TCR_Pielou_evenness_index_boxplot.svg' ,width = 6, height = 6)

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
ggsave('Th2_subsets_TCR_div.svg' ,width = 8, height = 8)

get_pairs <- function(char){
    char <- unique(char)
    result_x <- lapply(char, function(x){
        if(which(char == x) != length(char)){
            lapply(char[(which(char == x)+ 1):length(char)], function(y){
                return(c(x, y))
            })            
        }
    })
    result_x <- unlist(result_x, recursive = F)
    result_x <- result_x[!sapply(result_x, is.null)]
    return(result_x)
}

names(tcr_raw$data)

total_TRSS <- lapply(1:length(tcr_raw$data), function(i){
    nt_df <-   tcr_raw$data[[i]]
    #nt_df$identity[nt_df$identity %in% c('CD8+RTE','CD8+TN')] <- 'CD8+TN'
    
    #nt_df$tissue_identity <- paste0(nt_df$tissue, '_',nt_df$identity)
    pairs <- get_pairs(unique(nt_df$identity))
    pairs <- c(pairs, lapply(pairs, function(x) rev(x)))
    clu_val <- sapply(pairs, function(pair){
        clu_1_df <- nt_df[nt_df$identity %in% pair[1], ]
        clu_2_df <- nt_df[nt_df$identity %in% pair[2], ]
        clu_1_df_2 <- clu_1_df[clu_1_df$raw_clonotype_id %in% clu_2_df$raw_clonotype_id,]
        clu_2_df_2 <- clu_2_df[clu_2_df$raw_clonotype_id %in% clu_1_df$raw_clonotype_id,]  
        clu_1_stats <- table(clu_1_df_2$raw_clonotype_id)
        clu_2_stats <- table(clu_2_df_2$raw_clonotype_id)    
        val <- sapply(names(clu_1_stats), function(raw_clonotype_id){
            sqrt(clu_1_stats[raw_clonotype_id]/nrow(clu_1_df) * clu_2_stats[raw_clonotype_id]/nrow(clu_2_df))
        })
        val <- sum(unlist(val))              
    })    
    cluster1 <- sapply(pairs, function(x) x[1])
    cluster2 <- sapply(pairs, function(x) x[2])                   
    df <- data.frame(cluster1 = cluster1, cluster2 = cluster2, TRSS = clu_val, 
                     sample = names(tcr_raw$data)[i])        
      

    return(df)
})
total_TRSS <- do.call(rbind, total_TRSS)


total_TRSS$class <- paste0(total_TRSS$cluster1, '~', total_TRSS$cluster2)
TRSS_sum <- ddply(total_TRSS, `.variables` = c('class'), summarise, TRSS = mean(TRSS))
TRSS_sum$cluster1 <- sapply(strsplit(TRSS_sum$class, split = '~'), function(x){x[1]})
TRSS_sum$cluster2 <- sapply(strsplit(TRSS_sum$class, split = '~'), function(x){x[2]})

lbls_func <- function(x){
    x[x=='0'] <-  'C0: Stem-like'
    x[x=='1'] <- 'C1: Pathogenic'
    x[x=='2'] <- 'C2: Ikzf2+'
    x[x=='3'] <- 'C3: S100a4+'
    x[x=='4'] <- 'C4: Isg15+'
    x[x=='5'] <-  'C5: Mki67+'
    x
}



plot_df <- TRSS_sum#[TRSS_sum$cluster1 %in% cluster_used &  TRSS_sum$cluster2 %in% cluster_used, ]
plot_df <- rbind(plot_df, 
                 data.frame(class = paste0(0:5, '~',0:5), TRSS = NA, 
                            cluster1 = as.character(0:5),cluster2 = as.character(0:5) 
                           ))

min_val <- min(plot_df$TRSS)
max_val <- max(plot_df$TRSS)
#tmp_breaks <- round(seq(min_val, max_val, length.out = 5),2)

p <- ggplot(plot_df, aes(cluster1, cluster2))+
geom_point(shape = 21, mapping = aes(fill = TRSS, size = TRSS))+
geom_tile(fill =NA, color = 'black', size = 0.1)+
scale_fill_gradientn(colours = c("white","#F5E1DA","#EFCBBF","#A21F2A","#6F151C"),# breaks = tmp_breaks
                     #breaks = round(c(min_val, mean(c(min_val, max_val)), max_val), 2)
                    )+
scale_size_continuous(range = c(6,15), # breaks = tmp_breaks
                     #breaks = round(c(min_val, mean(c(min_val, max_val)), max_val), 2)
                     )+
scale_x_discrete(label = lbls_func)+
scale_y_discrete(limits = rev, label = lbls_func)+
theme(#plot.margin = margin(b = 0.7,l = 0.7, unit = 'cm'),
      axis.title.x = element_text(size=20),axis.ticks = element_blank(),strip.background = element_blank(),
      legend.position = 'right',
          axis.title.y = element_blank(),
         axis.text.x = element_text(size=20, hjust=1, angle = 45), 
      axis.text.y = element_text(size=20),
          axis.line  = element_blank(), 
         panel.border=element_rect(colour = "black", fill =NA),
          panel.background = element_rect(colour = NA, fill ='white'),
          legend.key.width = unit(0, "lines"),
       legend.key.height = unit(0, "lines"), legend.title = element_text(size=20),
          plot.title = element_text(size=20, hjust=0.5),
          legend.text = element_text(size=20),
      strip.text.x = element_text(size=20))+
labs(x = NULL, y = NULL)+
# geom_tile(data = plot_df[plot_df$class %in% c('0~0','1~1','2~2','3~3','4~4','5~5'), ], color = 'black',
#           fill = 'grey')+
guides(fill = guide_legend(shape = 19, color = 'black', ncol = 1))

options(repr.plot.width=9, repr.plot.height=7)
p
ggsave('TRSS.svg', width = 9, height = 7)
