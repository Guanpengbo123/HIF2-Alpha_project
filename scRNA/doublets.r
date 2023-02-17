library(Seurat, lib.loc = .libPaths()[2])
setwd('/data1/02.private/dengyj/analysis/Hif2a/doublets')

sub_dir <- '/filtered_feature_bc_matrix/'

dir1 <- '/data1/02.private/dengyj/analysis/Hif2a/matrix/'
sample_dir_1 = dir(dir1)
sample_dir_1 = sample_dir_1[grep('10X5', sample_dir_1)]
sample_1 <- sample_dir_1
sample_dir_1 <- paste0(dir1, sample_dir_1, sub_dir)

dir3 ='/data3/Hif2a/ABFC20200274-161728/ABFC20200274_10X-scRNA_result/2.Cellranger/'
sample_dir_3 = dir(dir3)
sample_dir_3 = sample_dir_3[grep('10X5', sample_dir_3)]
sample_3 <- sample_dir_3
sample_dir_3 <- paste0(dir3, sample_dir_3, sub_dir)

dir4 = '/data3/Hif2a/ABFC20200274--303137-10X大测分析结果/result/'
sample_dir_4 = dir(dir4)
sample_dir_4 = sample_dir_4[grep('10X5', sample_dir_4)]
sample_4 <- sample_dir_4
sample_dir_4 <- paste0(dir4, sample_dir_4, sub_dir)

sample_dir <- c(sample_dir_1,sample_dir_3,sample_dir_4)
sample <- c(sample_1, sample_3, sample_4)


for(i in 1:length(sample_dir)){
  doublets_list <- read.table(paste0('/data1/02.private/dengyj/analysis/Hif2a/doublets/', sample[i], '_doublets.txt'))
  iter_matrix<- Read10X(sample_dir[i])[, doublets_list$V1 ==0]
  iter_object <- CreateSeuratObject(iter_matrix, project = sample[i], min.cells = 3, min.features = 200)
  iter_object  <- RenameCells(iter_object, new.names = paste(colnames(iter_object), sample[i], sep = '_'))
  iter_object$batch <- sample[i]
  assign(sample[i], iter_object)
  print(paste0('Finished--', sample[i]))
}

hif_sample_list <- mget(sample)
hif_step_1 <- merge(hif_sample_list[[1]], hif_sample_list[2:length(hif_sample_list)])
save(hif_step_1, file = 'hif_step_1.rds')

