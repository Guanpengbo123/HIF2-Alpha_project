#####寻找共定位基因###
###linux

bedtools intersect -a GATA3.rmbl.bed -b HIF2a.rmbl.bed -wa -wb >m.bed
##R
common<-as.data.frame(fread("/data4/cut_TAG/HIF2_Gata3/SEACR_callpeak/Top_005_new/annofile/m.bed"))
common_FT<-common[common$V1%in%c(1:22,"X","Y"),]
for (i in (1:nrow(common_FT))){
    a<-common_FT$V2[i]
    b<-common_FT$V3[i]
    c<-common_FT$V8[i]
    d<-common_FT$V9[i]
    common_FT$chr[i]<-paste0("chr",common_FT$V1[i])
    common_FT$strat[i]<-max(a,c)
    common_FT$end[i]<-min(b,d)
    
}
write.table(common_FT[,c("chr","strat","end")],"./common_peaks_ft.bed",sep = "\t",col.names = F,row.names = F)
######Linux
bedtools intersect -a GATA3.rmbl.bed -b common_peaks_ft.bed -v > Only_GATA3_peaks.bed

bedtools intersect -a ./HIF2a.rmbl.bed -b common_peaks_ft.bed -v > Only_HIF2a_peaks.bed


#####基因注释
#########R
library("ChIPseeker")
library("org.Hs.eg.db")
library("clusterProfiler")
library("GenomicFeatures")
library(corrplot)
library(enrichplot)
library(ReactomePA)
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(dplyr)
libary(data.table)
GATA<-readPeakFile("/data4/cut_TAG/HIF2_Gata3/SEACR_callpeak/Top_005_new/GATA3_ft.bed")
HIF<-readPeakFile("/data4/cut_TAG/HIF2_Gata3/SEACR_callpeak/Top_005_new/HIF2_ft.bed")
GATA_Anno <- annotatePeak(GATA, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
HIF_Anno <- annotatePeak(HIF, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
GATA_adata<-as.data.frame(GATA_Anno)
HIF_adata<-as.data.frame(HIF_Anno)
GATA_adata$anno<-""
for (i in 1:length(rownames(GATA_adata))){
    a<-GATA_adata$annotation[i]
    b<-unlist(strsplit(a,''))
    n<-length(which(b=="("))
    if (n==0){
        GATA_adata$anno[i]<-a
    }else{
        n1<-which(b=="(")
        c<-b[1:n1-1]
        GATA_adata$anno[i]<-c(paste(c,collapse=''))
    }
}
HIF_adata$anno<-""
for (i in 1:length(rownames(HIF_adata))){
    a<-HIF_adata$annotation[i]
    b<-unlist(strsplit(a,''))
    n<-length(which(b=="("))
    if (n==0){
        HIF_adata$anno[i]<-a
    }else{
        n1<-which(b=="(")
        c<-b[1:n1-1]
        HIF_adata$anno[i]<-c(paste(c,collapse=''))
    }
}
GH_common_gene<-intersect(unique(GATA_adata$SYMBOL),unique(HIF_adata$SYMBOL))
only_GATA3_gene<-setdiff(GATA_adata$SYMBOL,GH_common_gene)
only_HIF_gene<-setdiff(HIF_adata$SYMBOL,GH_common_gene)

wb <- createWorkbook()
for (i in sample){ 
    genelist<-get(i)
    doc = read_docx()
    ego <- enrichGO(genelist, 
                    keyType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db, 
                    ont = "all", 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    minGSSize=5,
                    maxGSSize=500,
                    readable = TRUE)

    name1<-paste0(i,"_GO")
    
#dotplot(ego,showCategory=50)
#write.table(as.data.frame(ego),file=paste0("/data4/cut_TAG/HIF2_Gata3/annotext/pathway319/",i,"_go.xlsx"),sep="\t",col.names=T,row.names =F,quote = F)
    doc = read_docx()
    kegg <- enrichKEGG(genelist,
                 organism ="hsa",
                 pvalueCutoff = 1,
                 qvalueCutoff = 1,
                minGSSize = 5,
                  maxGSSize = 500,
                 use_internal_data =FALSE)
    KEGG<- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    name2<-paste0(i,"_KEGG")
#dotplot(kegg, showCategory = 20)
#write.table(as.data.frame(KEGG),file=paste0("/data4/cut_TAG/HIF2_Gata3/annotext/pathway319/",i,"_kegg.xlsx"),sep="\t",col.names=T,row.names =F,quote = F)
    doc = read_docx()
    Reactome<-enrichPathway(genelist,
                        readable=T,
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                            minGSSize = 5,
  maxGSSize = 500,
                       )
    name3<-paste0(i,"_Rtm")
#dotplot(Reactome,showCategory=20)
#write.table(as.data.frame(Reactome),file=paste0("/data4/cut_TAG/HIF2_Gata3/annotext/pathway319/",i,"_Rtm.xlsx"),sep="\t",col.names=T,row.names =F,quote = F)
    
    all_data<-rbind(as.data.frame(ego)[1:100,c("Description","pvalue","GeneRatio","BgRatio","Count")],
                    as.data.frame(KEGG)[1:100,c("Description","pvalue","GeneRatio","BgRatio","Count")],
                    as.data.frame(Reactome)[1:100,c("Description","pvalue","GeneRatio","BgRatio","Count")]
                   )
    all_data_ft<-all_data[all_data$pvalue<=0.05&all_data$Count>=10,]
    all_data_choose<-all_data_ft[,c("Description","pvalue","GeneRatio","BgRatio")]

    addWorksheet(wb, i)
    writeData(wb, sheet = i, 
          x = all_data_choose)
}
saveWorkbook(wb, file = "HIF2_GATA3_pathway_gene_analysis_ORA_0217_new.xlsx", overwrite = TRUE)


#####common_peak###

common_peak<-readPeakFile("/data4/cut_TAG/HIF2_Gata3/SEACR_callpeak/Top_005_new/annofile/common_peaks_ft.bed")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
common_peak_Anno <- annotatePeak(common_peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
common_adata<-as.data.frame(common_peak_Anno)
common_adata$anno<-""
for (i in 1:length(rownames(common_adata))){
    a<-common_adata$annotation[i]
    b<-unlist(strsplit(a,''))
    n<-length(which(b=="("))
    if (n==0){
        common_adata$anno[i]<-a
    }else{
        n1<-which(b=="(")
        c<-b[1:n1-1]
        common_adata$anno[i]<-c(paste(c,collapse=''))
    }
}

common_gene<-unique(common_adata$geneId)
common_pgene<-unique(common_adata$geneId[common_adata$anno=="Promoter "])


ff<-Th2_C01[rownames(Th2_C01)%in%commongene_Mmgene,]
ff
#ff_HSC_TSS<-ft[rownames(ft)%in%genelist_GSEA_TSS,]
#ff_HSC_TSS


get_gsea_file <- function(
    object, 
    idents_vars = NULL, 
    cls_name = 'cls.cls', 
    exp_name = 'matrix.gct',
    gseapy = T, 
    my_sep = '\t', 
    add_prefix=T, 
    prefix = 'C_'
){  
    write_cls <- function(x, object, cls_name){
      sink(file = cls_name)
      cat(ncol(object), my_sep, length(unique(x)), my_sep, 1, '\n', sep = '')
      cat('#', my_sep, paste(unique(x), collapse = my_sep), '\n', sep = '')
      cat(paste(x, collapse = my_sep), '\n', sep = '')
      sink()
    }
    if(!is.null(idents_vars)){
      Idents(object) <- object[[idents_vars]]
    }    
    clu_num <- length(levels(object)) 
    list <- NULL
    clu_name <- NULL
    lbl <- NULL
    for( i in 1:clu_num) {
      list2 <- WhichCells(object, idents = levels(object)[i])
      list <- if(is.null(list)) {list2} else {c(list, list2)}
      lbl2 <- if(add_prefix){
        rep(paste0('C_', levels(object)[i]), length(list2))
      } else {
        rep(levels(object)[i], length(list2))
      }
      lbl <- if(is.null(lbl)) {lbl2} else {c(lbl, lbl2)}
    }
    matrix <- as.data.frame(object$RNA@data[, list])
    matrix <- cbind(data.frame(Name = rownames(matrix), Description = rep("", nrow(matrix))), matrix)
    if(gseapy){
      cls <- matrix(lbl, nrow=1, ncol=length(lbl))
      write.table(cls, file = cls_name, quote = F, col.names = F, row.names = F, sep = my_sep)
    }else {
       write_cls(x = lbl, object = object, cls_name = cls_name)
    }    
    write.table(matrix, quote = F, sep = my_sep, row.names = F, file = exp_name)
}


get_gsea_file(ff, idents_vars = 'finally_class', exp_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/scRNA/HIF2/GSEA/input/GH_common_gene_SEACR005_C01_cluster_GSEA_matrix.txt', 
              cls_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/scRNA/HIF2/GSEA/input/GH_common_gene_SEACR005_C01_cluster_GSEA_matrix.cls',gseapy = F, add_prefix = T)

ff<-Th2_C01[rownames(Th2_C01)%in%onlygene_GATA3_Mmgene,]
ff
#ff_HSC_TSS<-ft[rownames(ft)%in%genelist_GSEA_TSS,]
#ff_HSC_TSS


get_gsea_file <- function(
    object, 
    idents_vars = NULL, 
    cls_name = 'cls.cls', 
    exp_name = 'matrix.gct',
    gseapy = T, 
    my_sep = '\t', 
    add_prefix=T, 
    prefix = 'C_'
){  
    write_cls <- function(x, object, cls_name){
      sink(file = cls_name)
      cat(ncol(object), my_sep, length(unique(x)), my_sep, 1, '\n', sep = '')
      cat('#', my_sep, paste(unique(x), collapse = my_sep), '\n', sep = '')
      cat(paste(x, collapse = my_sep), '\n', sep = '')
      sink()
    }
    if(!is.null(idents_vars)){
      Idents(object) <- object[[idents_vars]]
    }    
    clu_num <- length(levels(object)) 
    list <- NULL
    clu_name <- NULL
    lbl <- NULL
    for( i in 1:clu_num) {
      list2 <- WhichCells(object, idents = levels(object)[i])
      list <- if(is.null(list)) {list2} else {c(list, list2)}
      lbl2 <- if(add_prefix){
        rep(paste0('C_', levels(object)[i]), length(list2))
      } else {
        rep(levels(object)[i], length(list2))
      }
      lbl <- if(is.null(lbl)) {lbl2} else {c(lbl, lbl2)}
    }
    matrix <- as.data.frame(object$RNA@data[, list])
    matrix <- cbind(data.frame(Name = rownames(matrix), Description = rep("", nrow(matrix))), matrix)
    if(gseapy){
      cls <- matrix(lbl, nrow=1, ncol=length(lbl))
      write.table(cls, file = cls_name, quote = F, col.names = F, row.names = F, sep = my_sep)
    }else {
       write_cls(x = lbl, object = object, cls_name = cls_name)
    }    
    write.table(matrix, quote = F, sep = my_sep, row.names = F, file = exp_name)
}


get_gsea_file(ff, idents_vars = 'finally_class', exp_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/scRNA/HIF2/GSEA/input/only_GATA3_gene_SEACR005_C01_cluster_GSEA_matrix.txt', 
              cls_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/scRNA/HIF2/GSEA/input/only_GATA3_gene_SEACR005_C01_cluster_GSEA_matrix.cls',gseapy = F, add_prefix = T)

ff<-Th2_C01[rownames(Th2_C01)%in%onlygene_HIF_Mmgene,]
ff
#ff_HSC_TSS<-ft[rownames(ft)%in%genelist_GSEA_TSS,]
#ff_HSC_TSS


get_gsea_file <- function(
    object, 
    idents_vars = NULL, 
    cls_name = 'cls.cls', 
    exp_name = 'matrix.gct',
    gseapy = T, 
    my_sep = '\t', 
    add_prefix=T, 
    prefix = 'C_'
){  
    write_cls <- function(x, object, cls_name){
      sink(file = cls_name)
      cat(ncol(object), my_sep, length(unique(x)), my_sep, 1, '\n', sep = '')
      cat('#', my_sep, paste(unique(x), collapse = my_sep), '\n', sep = '')
      cat(paste(x, collapse = my_sep), '\n', sep = '')
      sink()
    }
    if(!is.null(idents_vars)){
      Idents(object) <- object[[idents_vars]]
    }    
    clu_num <- length(levels(object)) 
    list <- NULL
    clu_name <- NULL
    lbl <- NULL
    for( i in 1:clu_num) {
      list2 <- WhichCells(object, idents = levels(object)[i])
      list <- if(is.null(list)) {list2} else {c(list, list2)}
      lbl2 <- if(add_prefix){
        rep(paste0('C_', levels(object)[i]), length(list2))
      } else {
        rep(levels(object)[i], length(list2))
      }
      lbl <- if(is.null(lbl)) {lbl2} else {c(lbl, lbl2)}
    }
    matrix <- as.data.frame(object$RNA@data[, list])
    matrix <- cbind(data.frame(Name = rownames(matrix), Description = rep("", nrow(matrix))), matrix)
    if(gseapy){
      cls <- matrix(lbl, nrow=1, ncol=length(lbl))
      write.table(cls, file = cls_name, quote = F, col.names = F, row.names = F, sep = my_sep)
    }else {
       write_cls(x = lbl, object = object, cls_name = cls_name)
    }    
    write.table(matrix, quote = F, sep = my_sep, row.names = F, file = exp_name)
}


get_gsea_file(ff, idents_vars = 'finally_class', exp_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/scRNA/HIF2/GSEA/input/only_HIF_gene_SEACR005_C01_cluster_GSEA_matrix.txt', 
              cls_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/scRNA/HIF2/GSEA/input/only_HIF_gene_SEACR005_C01_cluster_GSEA_matrix.cls',gseapy = F, add_prefix = T)
