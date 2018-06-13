setClass(
  Class = "KEGG.annotater", 
  representation = representation(de.genes = 'data.frame',
                                  all.genes = 'data.frame',
                                  species = 'character',
                                  kegg.annotation = 'data.frame',
                                  de.kegg.genes = 'data.frame',
                                  de.kegg.genes.filtered = 'data.frame',
                                  all.kegg.genes = 'data.frame',
                                  all.kegg.genes.filtered = 'data.frame')
  
)

setGeneric("get.kegg.annotation", function(.Object) 
  standardGeneric("get.kegg.annotation") )

setGeneric("couple.kegg.annotation", function(.Object) 
  standardGeneric("couple.kegg.annotation") )

setGeneric("filter.kegg.pathways", function(.Object, priority.interest) 
  standardGeneric("filter.kegg.pathways") )

setGeneric("plot.fold.changes", function(.Object, type) 
  standardGeneric("plot.fold.changes") )

setGeneric("map.kegg.pathway", function(.Object, path.ids) 
  standardGeneric("map.kegg.pathway") )

setGeneric("remove.junk", function(.Object, path.id) 
  standardGeneric("remove.junk") )


#' @description This forms the init of the class
#' @param all.genes, a data frame with all the genes obtained from the SeqDat class
#' @param de.genes, a data frame with the DE genes obtained from the SeqDat class
setMethod("initialize", "KEGG.annotater", function(.Object, all.genes, de.genes, species){
  .Object@all.genes <- all.genes
  .Object@de.genes <- de.genes
  .Object@species <- species
  .Object
})

#' @description This method queries the KEGG database for the species name and 
#' retrieves the corresponding pathway and gene information
setMethod("get.kegg.annotation", signature("KEGG.annotater"), function(.Object) {
  # - Retrieving data from KEGG and make annotation table hereof. Note that
  # - we could have used the KEGGREST package. However this package
  # - returns a character vector insted of a data frame, which would require 
  # - an extra processing step. 
  pathway.id.name <- read.table(paste0("http://rest.kegg.jp/list/pathway/", .Object@species), quote="", sep="\t")
  gene.pathway.id <- read.table(paste0("http://rest.kegg.jp/link/pathway/",.Object@species), quote="", sep="\t")
  gene.info <- left_join(gene.pathway.id, pathway.id.name, by = c('V2'='V1'))
  
  # - Formatting makes downstream analysis much easier
  colnames(gene.info) <- c("ORF",'pathway.id','pathway.name')
  gene.info$ORF <- gsub(paste0(.Object@species,':'),'',gene.info$ORF)
  gene.info$pathway.id <- gsub('path:','',gene.info$pathway.id)
  gene.info$pathway.name <- gsub(' -.*',"", gene.info$pathway.name)

  # - Saving the results in the Object
  .Object@kegg.annotation <- .Object@de.genes %>% left_join(gene.info)
  .Object
})


#' @description This function  couples the KEGG annotation to DE genes and to
#' all genes. 
setMethod("couple.kegg.annotation", signature("KEGG.annotater"), function(.Object) {
  .Object@de.kegg.genes <- .Object@de.genes %>% left_join(.Object@kegg.annotation)
  .Object@all.kegg.genes <- .Object@all.genes %>% left_join(.Object@kegg.annotation)
  .Object
})


#' @description This funtion filters the DE genes and all genes for pathways of the users interest
#' Note however that the old data is still saved. 
setMethod("filter.kegg.pathways", signature("KEGG.annotater"), function(.Object, priority.interest) {
  .Object@de.kegg.genes.filtered <- .Object@de.kegg.genes %>%
    filter(pathway.name  %in% priority.interest) %>%
    arrange(-logFC)
  .Object@all.kegg.genes.filtered <- .Object@all.kegg.genes %>%
    filter(pathway.name  %in% priority.interest) %>%
    arrange(-logFC)
  .Object
})


#' @description This function can be used to plot the spread of FC in the pathways of 
#' the users interest for de DE genes. 
#' @param type typ 1 will produce a graph for each pathway seperately, whereas type 2
#' will show the data for all pathways together.
setMethod("plot.fold.changes", signature("KEGG.annotater"), function(.Object, type) {
  if (type == 1) {
    # - A quick look at the genes and pathways
    ggplot(.Object@de.kegg.genes.filtered, aes(x = gsub(' ','\n',pathway.name), y=logFC, col = pathway.name)) + 
      geom_point(size = 2.5, show.legend = F) +
      xlab('Pathway') +
      ylab('LogFC') +
      ggtitle('Global pathways changes') +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_hline(yintercept = 0, col = 'red', size = 0.9) 
  } else {
    # - Alternatively, we could plot these genes ordered by FC and color according to the pathway
    ggplot(.Object@de.kegg.genes.filtered, aes(x = seq(1,nrow(.Object@de.kegg.genes.filtered)), y = logFC, col = pathway.name)) +
      geom_point(size = 2.5)
  }
})

#' @description This function  couples the KEGG annotation to DE genes and to
#' all genes. 
setMethod("map.kegg.pathway", signature("KEGG.annotater"), function(.Object, path.ids) {
  # Formatting the data for KEGG analysis
  kegg.format <- wcfs1@edger.de %>% dplyr::select(logFC)
  
  # Getting names for the provied path ids
  path.names <- .Object@kegg.annotation %>% 
    filter(pathway.id %in% path.ids) %>%
    .$pathway.name %>% unique()
  
  print(path.names)
  print(path.ids)
  
  # Helper function to map genes 
  map <- function(path.id, out.suffix) {
    pathview(gene.data = kegg.format, 
             gene.idtype = "KEGG" , 
             pathway.id = path.id , 
             species = SPECIES,
             map.symbol = TRUE,        # show gene names when possible
             same.layer = FALSE)       # extra layer for the gene names
  }
  mapply(map, path.ids, path.names)
})


#' @description This function  couples the KEGG annotation to DE genes and to
#' all genes. 
#' Pathview uses two files to create the output image:
#' PNG image of the KEGG pathay
#' XML file of the KEGG pathway
#' Those aren't deleted automatically, therefore we 
#' wrote this function
setMethod("remove.junk", signature("KEGG.annotater"), function(.Object, path.id) {
  remove.junk <- function(path.id) {
    png.base = paste0(path.id, '.png')
    xml.file = paste0(path.id, '.xml')
    if (file.exists(png.base)) file.remove(png.base)
    if (file.exists(xml.file)) file.remove(xml.file)
  }
})




















