#'###############################################################################################
#' This class RNA.seq.analyzer will only require a count matrix and the gene annotation
#' of the genes in the count matrix. The count matrix should be named using 
#' x.1 x.2 vs y.1 y.2 as the names will be used to set up the contrast. Notably, 
#' the contrast can easily be changed by simply re-ordering the columns in the count
#' matrix. Hereafter, the user can run either limma, DESeq2 or EdgeR using the defaults
#' (set up by us). Of course the parameters can be changed in this code. This class
#' also allows to construct a venn diagram showing the differences between these methods.
#' Lastly, this class saves the DE genes for each package using the same format, facilitating
#' comparison and allowing fair p-value comparison (as e.g. the same correction methods 
#' will be applied). The data from the venn diagram can be used to obtain all genes identiefied
#' as DE in all packages. 
#' 
#' @details: Last update 4-07-18
#' @author Rick Beeloo, Koen v.d. Heide and Thomas Reinders
#'###############################################################################################
setClass(
  Class = "RNA.seq.analyzer", 
  representation = representation(count.matrix = "matrix",
                                  sample.annotation = "data.frame",
                                  deseq2.data = "DESeqDataSet",
                                  deseq2.result = "data.frame",
                                  deseq2.de = "data.frame",
                                  edger.data = "DGEList",
                                  edger.result = "data.frame",
                                  edger.de = "data.frame",
                                  limma.data = "EList",
                                  limma.result = "data.frame",
                                  gene.annotation = "data.frame",
                                  limma.de = "data.frame",
                                  venn.data = 'list'))

#' @description This forms the init of the class, and automatically builds the sample
#' annotation (later used for the contrast and design) based on the column names
#' in the count matrix
#' @param count.matrix, a matrix containing the RNA-seq counts in which each column
#' indicated a sample and each row a gene
#' @param gene.annotation, a data frame in which each row is a gene, and the first 
#' row should be named (ORF) all other metadata columns are not directly used and therefore
#' not restricted to naming.  
setMethod("initialize", "RNA.seq.analyzer", function(.Object, count.matrix, gene.annotation){
  .Object@count.matrix <- count.matrix
  .Object@gene.annotation <- gene.annotation
  .Object@sample.annotation = data.frame(condition = gsub("\\.1|.2", "", colnames(count.matrix)))
  .Object
})


setGeneric("run.limma", function(.Object) 
  standardGeneric("run.limma") )

setGeneric("run.edger", function(.Object) 
  standardGeneric("run.edger") )

setGeneric("run.deseq2", function(.Object) 
  standardGeneric("run.deseq2") )

setGeneric("filter.low.counts", function(.Object, cpm.cut.off) 
  standardGeneric("filter.low.counts") )

setGeneric("determine.significance", function(.Object, p.val) 
  standardGeneric("determine.significance") )

setGeneric("draw.venn.diagram", function(.Object) 
  standardGeneric("draw.venn.diagram") )

setGeneric("get.annotation", function(.Object, wanted.ids) 
  standardGeneric("get.annotation") )

setGeneric("run.all", function(.Object) 
  standardGeneric("run.all") )

#' @description This method filters out genes with low counts
#' @param wanted.ids, a list of gene ids for which the annotation should 
#' be retrieved.
setMethod("filter.low.counts", signature("RNA.seq.analyzer"), function(.Object, cpm.cut.off) {
  count.data <- .Object@count.matrix
  keep <- rowSums(cpm(count.data)>1) >= 2
  .Object@count.matrix <- count.data[keep,]
  .Object
})

#' @description This method can be used to run limma in default settings
#' The original results will be saved in Object@limma.data (EList object), 
#' whereas the standarized results will be saved in Object@limma.result
setMethod("run.limma", signature("RNA.seq.analyzer"), function(.Object) {
  require(limma)
  require(edgeR)
  groups = factor(.Object@sample.annotation$condition)
  nf <- edgeR::calcNormFactors(.Object@count.matrix, method = 'TMM')
  design = model.matrix(~groups)
  voom.data <- limma::voom(.Object@count.matrix, 
                           design,
                           lib.size = colSums(.Object@count.matrix) * nf)
  voom.data$genes <- rownames(.Object@count.matrix)
  voom.fitlimma <- limma::lmFit(voom.data, design = design)
  voom.fitbayes <- limma::eBayes(voom.fitlimma)
  voom.pvalues <- voom.fitbayes$p.value[, 2]
  voom.adjpvalues <- p.adjust(voom.pvalues, method = 'fdr')
  voom.logFC <- voom.fitbayes$coefficients[, 2]
  result.table <- data.frame('pvalue' = voom.pvalues, 'adjpvalue' = voom.adjpvalues, 'logFC' = voom.logFC)
  rownames(result.table) <- rownames(.Object@count.matrix)
  print(paste0("contrast used: ", unique(groups)[1], ' vs ', unique(groups)[2]))
  .Object@limma.result <- result.table
  .Object@limma.data <- voom.data
  .Object
})

#' @description This method can be used to run edgeR in default settings
#' The original results will be saved in Object@edger.data (DGEList object), 
#' whereas the standarized results will be saved in Object@edger.result 
setMethod("run.edger", signature("RNA.seq.analyzer"), function(.Object) {
  require(edgeR)
  groups = factor(.Object@sample.annotation$condition)
  edgeR.dgelist <- edgeR::DGEList(counts = .Object@count.matrix, group = groups)
  edgeR.dgelist <- edgeR::calcNormFactors(edgeR.dgelist, method = 'TMM')
  design <- model.matrix(~0+group, data = edgeR.dgelist$samples)
  colnames(design) <- levels(edgeR.dgelist$samples$group)
  edgeR.dgelist <- edgeR::estimateGLMCommonDisp(edgeR.dgelist, design = design)
  edgeR.dgelist <- edgeR::estimateGLMTrendedDisp(edgeR.dgelist, design = design, method = 'power')
  edgeR.dgelist <- edgeR::estimateGLMTagwiseDisp(edgeR.dgelist, design = design)
  edgeR.fit <- edgeR::glmFit(edgeR.dgelist , design)
  def.con <<- paste0(unique(.Object@sample.annotation$condition)[2],
                '-', unique(.Object@sample.annotation$condition)[1])
  mc <- limma::makeContrasts(def.con , levels = design)
  edgeR.glrt <- edgeR::glmLRT(edgeR.fit, contrast = mc)
  edgeR.dgelist$glmrt <- edgeR::glmLRT(edgeR.fit, contrast = mc)
  edgeR.pvalues <- edgeR.glrt$table$PValue
  edgeR.adjpvalues <- p.adjust(edgeR.pvalues, method = 'fdr')
  edgeR.logFC <- edgeR.glrt$table$logFC
  result.table <- data.frame('pvalue' = edgeR.pvalues, 'adjpvalue' = edgeR.adjpvalues, 'logFC' = edgeR.logFC)
  rownames(result.table) <- rownames(.Object@count.matrix)
  print(paste0("contrast used: ", unique(groups)[1], ' vs ', unique(groups)[2]))
  .Object@edger.result <- result.table
  .Object@edger.data <- edgeR.dgelist
  .Object
})

#' @description This method can be used to run DESeq2 in default settings
#' The original results will be saved in Object@deseq2.data (DERNA.seq.analyzeraSet), 
#' whereas the standarized results will be saved in Object@deseq2.result
setMethod("run.deseq2", signature("RNA.seq.analyzer"), function(.Object) {
  require(DESeq2)
  groups <- factor(.Object@sample.annotation$condition)
  DESeq2.ds <-DESeq2::DESeqDataSetFromMatrix(countData = .Object@count.matrix, 
                                              colData = data.frame(condition = groups), 
                                              design = ~ condition)
  DESeq2.ds <- DESeq2::DESeq(DESeq2.ds, quiet = TRUE)
  DESeq2.results <- DESeq2::results(DESeq2.ds, pAdjustMethod = 'fdr')
  DESeq2.pvalues <- DESeq2.results$pvalue
  DESeq2.adjpvalues <- DESeq2.results$padj
  DESeq2.logFC <- DESeq2.results$log2FoldChange
  result.table <- data.frame('pvalue' = DESeq2.pvalues, 'adjpvalue' = DESeq2.adjpvalues, 'logFC' = DESeq2.logFC)
  rownames(result.table) <- rownames(.Object@count.matrix)
  print(paste0("contrast used: ", unique(groups)[1], ' vs ', unique(groups)[2]))
  .Object@deseq2.result <- result.table
  .Object@deseq2.data <- DESeq2.ds
  .Object
})

#' @description This method can be used to determine the DE genes using a p.val cut-off
#' The results will be saved in Object@limma.de, Object@edger.de, Object@
#' deseq2.de (depending on which of these was/were used)
#' @param p.val, the p-value cut-off that should be used to determine
#' significance. 
setMethod("determine.significance", signature("RNA.seq.analyzer"), function(.Object, p.val) {
  if (nrow(.Object@limma.result) != 0) {
    de<- .Object@limma.result[.Object@limma.result$adjpvalue < p.val,]
    print(paste0('limma: ', nrow(de), ' significant genes'))
    .Object@limma.de <- de
  }
  if (nrow(.Object@edger.result) != 0) {
    de <- .Object@edger.result[.Object@edger.result$adjpvalue < p.val,]
    print(paste0('edgeR: ', nrow(de), ' significant genes'))
    .Object@edger.de <- de
  }
  if (nrow(.Object@deseq2.result) != 0) {
    de <- .Object@deseq2.result[.Object@deseq2.result$adjpvalue < p.val,]
    print(paste0('DESeq2: ', nrow(de), ' significant genes'))
    .Object@deseq2.de <- de
  }
  .Object
})

#' @description This method can be used to draw a venn diagram of the DE genes
#' This rquires that all the packages are called (either seperatetely 
#' or all together using run.all) 
setMethod("draw.venn.diagram", signature("RNA.seq.analyzer"), function(.Object) {
  require(VennDiagram)
  if (nrow(.Object@limma.de) != 0 & 
      nrow(.Object@edger.de) != 0 &
      nrow(.Object@deseq2.de) != 0) {
    data <- list('limma' = rownames(.Object@limma.de), 
                 'edger' = rownames(.Object@edger.de), 
                 'deseq2' = rownames(.Object@deseq2.de))
    while (!is.null(dev.list()))  dev.off() # only turn it off when its on
    v <- venn.diagram(data, filename = NULL, 
                      alpha = 0.5, fill=c("darkmagenta", "darkblue",'red'))
    grid.draw(v)
    .Object@venn.data <- get.venn.dat(data)
  } else {
    print('Run the "determine.significance" code first')
  }
  .Object
})

#' @description This method can be used to retrieve te annotation for a list of gene ids,
#' @param wanted.ids, a list of gene ids for which the annotation should 
#' be retrieved.
setMethod("get.annotation", signature("RNA.seq.analyzer"), function(.Object, wanted.ids) {
  anno <- .Object@gene.annotation
  anno <- anno[anno$ORF %in% wanted.ids,]
  anno
})

#' @description This method can be used to run all DESeq2, edgeR and limma at once
setMethod("run.all", signature("RNA.seq.analyzer"), function(.Object) {
  .Object <- run.deseq2(.Object)
  .Object <- run.edger(.Object)
  .Object <- run.limma(.Object)
  .Object
})
















