###code used to detect Th2 differentation via Monocle3
library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)
library(openxlsx)
library(monocle3)

Hif2a_Th2_2 <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/basic_analysis//Hif2a_Th2_2.rds')

Idents(Hif2a_Th2_2) <- Hif2a_Th2_2$seurat_clusters

hif_step_2 <- readRDS('/data1/02.private/dengyj/analysis/Hif2a/basic_analysis//hif_step_2.rds')


###seurat_clusters,0=pathogenic Th2, 2 = Treg-like Th2, 3 = stem-like Th2
Th2_dev <- merge(Hif2a_Th2_2[, Hif2a_Th2_2$seurat_clusters %in% c('0', '2', '3')], 
                 hif_step_2[, Idents(hif_step_2) == 'T_naive' & hif_step_2$source %in% c('WT', 'KO')])



Th2_dev_cds <- 
new_cell_data_set(Th2_dev$RNA@counts, 
                  gene_metadata = data.frame(gene_short_name = rownames(Th2_dev$RNA@counts), 
                                            row.names = rownames(Th2_dev$RNA@counts)), 
                  cell_metadata = Th2_dev@meta.data)

pData(Th2_dev_cds)$idents <- Idents(Th2_dev)

Th2_dev_cds <- preprocess_cds(Th2_dev_cds, num_dim =30)
Th2_dev_cds <- align_cds(Th2_dev_cds, alignment_group = "batch",residual_model_formula_str = NULL)


Th2_dev_cds <- reduce_dimension(Th2_dev_cds, umap.min_dist =0.5)


options(repr.plot.width = 9, repr.plot.height = 9)
plot_cells(Th2_dev_cds, color_cells_by = "idents")


options(repr.plot.width = 9, repr.plot.height = 9)
###looking for a Monocle3 clustering result that was similar to Seurat
for(i in 1:50){
    Th2_dev_cds <- cluster_cells(Th2_dev_cds, resolution = 5e-4, random_seed = i)
    p <- plot_cells(Th2_dev_cds, label_groups_by_cluster=FALSE,  color_cells_by = "cluster", group_label_size=8, 
           show_trajectory_graph = F,
          graph_label_size = 4)+
    ggtitle(paste0('seed: ', i))
    print(p)
}

Th2_dev_cds <- cluster_cells(Th2_dev_cds, resolution = 5e-4, random_seed = 14)

options(repr.plot.width = 9, repr.plot.height = 9)
plot_cells(Th2_dev_cds, label_groups_by_cluster=FALSE,  color_cells_by = "cluster", group_label_size=8, 
           show_trajectory_graph = F,
          graph_label_size = 4)

Th2_dev_cds <- learn_graph(Th2_dev_cds)


options(repr.plot.width = 9, repr.plot.height = 9)
plot_cells(Th2_dev_cds, color_cells_by = "idents")

####to find resonable start point based on UMAP, this cell was a CD4 naive T-cell
which.min(reducedDims(Th2_dev_cds)[['UMAP']][, 2])

Th2_dev_cds <- order_cells(Th2_dev_cds, root_cells = 'TGGGCGTGTCCGTTAA-1_W4-10X5')

options(repr.plot.width = 9, repr.plot.height = 9)
plot_cells(Th2_dev_cds, label_groups_by_cluster=T,  color_cells_by = "pseudotime",
           group_label_size=8, cell_size = 0.5,
          graph_label_size = 5,label_cell_groups=F,
           label_leaves=TRUE,
           label_branch_points=TRUE)+
theme(legend.text = element_text(size = 15), legend.title = element_text(size = 20))+
labs(color = 'Identity')


saveRDS(Th2_dev, file ='Th2_dev.rds')

saveRDS(Th2_dev_cds, file ='Th2_dev_cds.rds')
