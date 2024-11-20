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



save(Hs_plot_df, Mm_plot_df, file = 'Epas1_Path2_cor.RData')

save(Hs_public_Th2_pools, Mm_public_Th2_pools, file = 'public_Th2_pools.RData')
