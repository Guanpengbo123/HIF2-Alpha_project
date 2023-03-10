library(Seurat)
library(future)
#library(GSVA)
library(openxlsx)
library(Matrix)
source('/data1/02.private/other/HSC/Cal_enrichment.R')
ssgsea <- function(
  exp, 
  genesets, 
  min_size = 10, 
  max_size = Inf, 
  alpha = 0.25, 
  exp_filter = F, 
  sd_cutoff = 0.05,
  normalization = F, 
  multiprocess = T)
{
  if(multiprocess){
    if(future::nbrOfWorkers() > 1){
      my_lapply <- future.apply::future_lapply
    }else{
      warning('Not enough workers for multi-processing, please run the future::plan function if you want to get more available workers.')
      my_lapply <- pbapply::pblapply
    }
  }else{
    my_lapply <- pbapply::pblapply
  }
  if(class(exp) != 'matrix'){
    if(class(exp) != 'dgCMatrix'){
      stop('The expression matrix you privided must be a matrix or dgCMatrix object')
    }
  }
  
  if(exp_filter){
    val <- apply(exp, 2, sd)
    exp <- exp[val >= sd_cutoff, ,drop = F]
  }
  
  filtered_genesets <- lapply(genesets, function(gs){
    which(rownames(exp) %in% gs)
  })
  names(filtered_genesets) <- names(genesets)
  gs_size <- sapply(filtered_genesets, length)
  filtered_genesets <- filtered_genesets[gs_size >= min_size & gs_size <= max_size]
  
  es = my_lapply(1:ncol(exp), function(i) {
    gene_rank <- as.integer(rank(exp[, i]))  ###gene_rank中，基因值越大，秩越大
    gene_rank_loc = order(gene_rank, decreasing = TRUE)###类似于GSEA的rank_list，表达越高的基因越靠前
    
    es_iter = sapply(filtered_genesets, function(geneset) {
      
      pos_loc = gene_rank_loc %in% geneset###
      pos_neg = !pos_loc
      
      rank_score_alpha  = (gene_rank[gene_rank_loc] * pos_loc) ^ alpha
      
      step_cdf_pos = cumsum(rank_score_alpha)    / sum(rank_score_alpha)
      step_cdf_neg = cumsum(pos_neg) / sum(pos_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
      return(sum(step_cdf_diff))
    })
    return(unlist(es_iter))
  })
  es <- do.call(cbind, es)
  if(normalization){
    es <- es/diff(range(es))
  }
  rownames(es) <- names(filtered_genesets)
  colnames(es) <- colnames(exp)
  return(es)
}

Th17 <- strsplit(readLines("/data1/02.private/guanpb/GSEA_GMT/GSE121599_Th17_genesets.gmt"), split = '\t')
geneset_names <- sapply(1:length(Th17), function(x) {Th17[[x]][1]})
Th17 <- sapply(1:length(Th17), function(x) {Th17[[x]][-c(1,2)]})
names(Th17) <- geneset_names

Treg_like <- strsplit(readLines("/data1/02.private/guanpb/GSEA_GMT/GSE121599_Th17_Treg_genesets.gmt"), split = '\t')
geneset_names <- sapply(1:length(Treg_like), function(x) {Treg_like[[x]][1]})
Treg_like <- sapply(1:length(Treg_like), function(x) {Treg_like[[x]][-c(1,2)]})
names(Treg_like) <- geneset_names

Hif2a_Th2<-readRDS("/data1/02.private/other/Hif2a_Th2_2.rds")
DefaultAssay(Hif2a_Th2) <- "RNA"
Hif2a_Th2 <- NormalizeData(Hif2a_Th2,verbose = T)
Hif2a_Th2 <- FindVariableFeatures(Hif2a_Th2)
matrix <- Hif2a_Th2$RNA@data[rownames(Hif2a_Th2$RNA@meta.features)[Hif2a_Th2$RNA@meta.features$vst.variance.standardized > 0.5], ]
plan('multisession', workers = 10)
Hif2a_Th2_es <- ssgsea(exp =matrix, genesets = Treg_like, min_size = 5, max_size = 5000, exp_filter = F, normalization = F)
plan('sequential')
