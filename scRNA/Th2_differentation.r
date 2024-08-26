library(data.table)
library(dplyr)
library(Seurat)

library(ggplot2)
library(openxlsx)
library(SCopeLoomR)
library(monocle)

library(monocle3)

library(harmony)

Hif2a_Th2_2 <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/basic_analysis//Hif2a_Th2_2.rds')

hif_step_2 <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/basic_analysis//hif_step_2.rds')

Idents(hif_step_2) <- hif_step_2$identity



options(repr.plot.width = 6, repr.plot.height = 6)
DimPlot(Hif2a_Th2_2, label =T, group.by = 'identity')



Idents(Hif2a_Th2_2) <- Hif2a_Th2_2$identity

Th2_dev <- merge(Hif2a_Th2_2[, !Hif2a_Th2_2$seurat_clusters %in% c('5')], 
                 hif_step_2[, hif_step_2$identity == 'T_naive' & hif_step_2$source %in% c('WT', 'KO')])



Th2_dev$idents <- Idents(Th2_dev)

Th2_dev$batch <- factor(Th2_dev$batch)
Th2_dev <- NormalizeData(Th2_dev)
Th2_dev <- FindVariableFeatures(Th2_dev, selection.method = 'vst')#, mean.cutoff = c(0.05, 8), 
                                #dispersion.cutoff = c(0.5, Inf))


Th2_dev <- ScaleData(Th2_dev, vars.to.regress = c('batch'))
Th2_dev <- RunPCA(Th2_dev)


options(repr.plot.width=5, repr.plot.height=5)
ElbowPlot(Th2_dev, ndims = 50, reduction = 'pca')

Th2_dev <- RunUMAP(Th2_dev, reduction = 'pca', dims = c(1:5,7:8), verbose = F)

options(repr.plot.width=8, repr.plot.height=6)
DimPlot(Th2_dev, group.by = 'idents')+
theme(panel.background = element_rect(color = NA, fill = 'white'), 
     legend.position = 'right',#c(0.15,0.75),
     panel.border = element_rect(fill = NA, colour = 'black'), 
      axis.line = element_blank(), axis.ticks = element_blank(), 
      axis.text.y = element_text(size = 15, face = 'italic'),
      axis.text.x = element_text(size = 15, angle = 45, hjust=1),
      legend.text = element_text(size = 15),
      legend.title =element_text(size = 15), 
      axis.title.x = element_blank(),
      axis.title.y = element_blank())+
guides(color = guide_legend(title = 'Identity', override.aes = list(size = 7)))


loom <- build_loom('Th2_dev.loom',dgem=Th2_dev$RNA@scale.data)
close_loom(loom)
write.csv(Embeddings(Th2_dev, reduction = 'pca')[,setdiff(1:50, c())], 
          quote = F, file = 'Th2_dev_pca.csv')
Th2_dev$idents <- as.character(Idents(Th2_dev))
write.csv(Th2_dev[[]], quote = F, file = 'Th2_dev_metadata.csv')

Th2_dev_cds <- 
new_cell_data_set(Th2_dev$RNA@counts, 
                  gene_metadata = data.frame(gene_short_name = rownames(Th2_dev$RNA@counts), 
                                            row.names = rownames(Th2_dev$RNA@counts)), 
                  cell_metadata = Th2_dev@meta.data)

pData(Th2_dev_cds)$idents <- Idents(Th2_dev)



reducedDims(Th2_dev_cds)$PCA <- Embeddings(Th2_dev, reduction = 'pca')
Th2_dev_cds@preprocess_aux$prop_var_expl <- Th2_dev@reductions$pca@stdev
reducedDims(Th2_dev_cds)$UMAP <- Embeddings(Th2_dev, reduction = 'umap')#Th2_dev_FDG_df[, c('x', 'y')]#
Th2_dev_cds@preprocess_aux$gene_loadings <- Th2_dev@reductions$pca@feature.loadings                                                   


options(repr.plot.width = 6, repr.plot.height = 6)
plot_cells(Th2_dev_cds, label_groups_by_cluster=FALSE,  cell_size = 0.5,
           color_cells_by = "identity", group_label_size=8, 
           show_trajectory_graph = F,
          graph_label_size = 4)


options(repr.plot.width = 12, repr.plot.height = 6)
for(i in 1:50){
    Th2_dev_cds <- cluster_cells(Th2_dev_cds, resolution = 3e-4, random_seed = i)
    p <- plot_cells(Th2_dev_cds, label_groups_by_cluster=FALSE,  color_cells_by = "cluster", group_label_size=8, 
           show_trajectory_graph = F,
          graph_label_size = 4)+
    ggtitle(paste0('seed: ', i))
    Th2_dev_cds <- learn_graph(Th2_dev_cds,close_loop = F)

    p1 <- plot_cells(Th2_dev_cds, label_groups_by_cluster=FALSE,  cell_size = 1,#label_groups_by_cluster = F,
           color_cells_by = "idents", group_label_size=8,   
                     label_leaves = F, label_roots = F,
                label_branch_points = F, 
          graph_label_size = 4)+
    ggtitle(paste0('seed: ', i))

    print(p+p1)

}

Th2_dev_cds <- cluster_cells(Th2_dev_cds, resolution = 3e-4, random_seed = 14)
Th2_dev_cds <- learn_graph(Th2_dev_cds,close_loop = F)

options(repr.plot.width = 12, repr.plot.height = 6)

p <- plot_cells(Th2_dev_cds, label_groups_by_cluster=FALSE,  color_cells_by = "cluster", group_label_size=8, 
       show_trajectory_graph = F,
      graph_label_size = 4)


p1 <- plot_cells(Th2_dev_cds, label_groups_by_cluster=FALSE,  cell_size = 1,#label_groups_by_cluster = F,
       color_cells_by = "idents", group_label_size=8,   
                 label_leaves = F, label_roots = F,
            label_branch_points = F, 
      graph_label_size = 4)

print(p+p1)


options(repr.plot.width = 12, repr.plot.height = 6)
p <- plot_cells(Th2_dev_cds, label_groups_by_cluster=FALSE,  color_cells_by = "cluster", group_label_size=8, 
       show_trajectory_graph = F,
      graph_label_size = 4)
p1 <- plot_cells(Th2_dev_cds, label_groups_by_cluster=FALSE,  cell_size = 1,#label_groups_by_cluster = F,
       color_cells_by = "idents", group_label_size=8,   
                 label_leaves = F, label_roots = F,
            label_branch_points = F, 
      graph_label_size = 4)
print(p+p1)


Th2_dev_cds <- order_cells(Th2_dev_cds, root_cells = 'GTTCTCGTCTCGTTTA-1_W3-10X5')
options(repr.plot.width = 6, repr.plot.height = 6)
plot_cells(Th2_dev_cds, label_groups_by_cluster=T,  color_cells_by = "pseudotime",
           group_label_size=8, cell_size = 0.5,
          graph_label_size = 5,label_cell_groups=F,
           label_leaves=TRUE,
           label_branch_points=TRUE)+
#guides(color = guide_legend(override.aes = list(size = 7)))+
theme(legend.text = element_text(size = 15), legend.title = element_text(size = 20))




cluster_colors <- 
c("Th2"="#E5D2DD","T_memory"="#53A85F","Th17"="#F1BB72","Treg"="#F3B1A0",
  "Alveolar_Mac"="#D6E7A3",
"IFN_T"="#57C3F3","Th1"="#476D87","Cd8_T"="#E95C59","Trans_mac"="#E59CC4","T_naive"="#AB3282",
"Neutrophil"="#23452F","Bcell"="#BD956A","Inflammatory_Mac"="#8C549C","Mki67_T"="#585658",
"Eosinophil"="#9FA3A8","Interstitial_mac"="#E0D4CA","Epi"="#5F3D69","cDC"="#C5DEBA",
  "γδT"="#58A4C3",
"Mast"="#E4C755","pDC"="#F7F398","17"="#AA9A59","8"="#E63863","other" = "gray80", 
 'Tcell'='#4284c9', 
 '1'='#d87aad', '0' = '#76b0ec', '2'='#2daaaa', '3'='#625D9E', '4' = '#68A180', 
  '5'='#3A6963')


options(repr.plot.width = 5, repr.plot.height = 5)
p <- plot_cells(Th2_dev_cds, label_groups_by_cluster=FALSE,  label_leaves = T, label_roots = T,
                label_branch_points = T,      
                cell_size = 1,#label_groups_by_cluster = F,
           color_cells_by = "identity", group_label_size=8,             
          graph_label_size = 4.5)+
scale_color_manual(values = cluster_colors, limits = c('T_naive', '0', '1', '2', '3', '4'))+
theme(panel.background = element_rect(color = NA, fill = 'white'), 
      legend.position = c(0.75,0.22),
     panel.border = element_rect(fill = NA, colour = 'black'), 
      axis.line = element_blank(), axis.ticks = element_blank(), 
      axis.text.y = element_text(size = 15, face = 'plain'),
      axis.text.x = element_text(size = 15, angle = 45, hjust=1),
      legend.text = element_text(size = 15),
      legend.title =element_text(size = 15), 
      axis.title.x = element_blank(),
      axis.title.y = element_blank())+
guides(color = guide_legend(title = 'Identity', ncol = 2, byrow = F,override.aes = list(size = 7)))
p

p1 <- p
p1$data <- p1$data[p1$data$source == 'WT', ]
options(repr.plot.width = 5, repr.plot.height = 5)
p1

p1 <- p
p1$data <- p1$data[p1$data$source == 'KO', ]
options(repr.plot.width = 5, repr.plot.height = 5)
p1

options(repr.plot.width = 5, repr.plot.height = 5)
plot_cells(Th2_dev_cds, label_groups_by_cluster=T,  color_cells_by = "pseudotime",label_leaves = T, 
           label_roots = T,
                label_branch_points = T,      
           group_label_size=8, cell_size = 0.5,
          graph_label_size = 4.5,label_cell_groups=F)+
#guides(color = guide_legend(override.aes = list(size = 7)))+
theme(panel.background = element_rect(color = NA, fill = 'white'), 
      legend.position = c(0.75,0.22),
     panel.border = element_rect(fill = NA, colour = 'black'), 
      axis.line = element_blank(), axis.ticks = element_blank(), 
      axis.text.y = element_text(size = 15, face = 'plain'),
      axis.text.x = element_text(size = 15, angle = 45, hjust=1),
      legend.text = element_text(size = 15),
      legend.title =element_text(size = 15), 
      axis.title.x = element_blank(),
      axis.title.y = element_blank())+
guides(color = guide_colorbar(title = 'Pseudotime', override.aes = list(size = 7)))
