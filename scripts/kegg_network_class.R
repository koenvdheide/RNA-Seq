setClass(
  Class = "KEGG.network", 
  representation = representation(de.genes = 'data.frame',
                                  nodes = 'data.frame',
                                  links = 'data.frame')
  
)

setGeneric("build.nodes", function(.Object) 
  standardGeneric("build.nodes") )

setGeneric("connect.nodes", function(.Object) 
  standardGeneric("connect.nodes") )

setGeneric("plot", function(.Object, name, direct.visible = FALSE) 
  standardGeneric("plot") )

#' @description This forms the init of the class
#' @params de.genes is a set of genes identified as differentially expressed with their
#' corresponding kegg annotation
setMethod("initialize", "KEGG.network", function(.Object, de.genes){
  .Object@de.genes <- de.genes
  .Object
})

#' @description This function builds a node data frame
#' all genes. 
setMethod("build.nodes", signature("KEGG.network"), function(.Object) {
  # - Building an ORF node set, we colored these nodes based on whether 
  # - the gene was identified as up or down regulated. Further 
  # - The annotation is added as it more intuitive to interpret gene
  # - names rather than the locus tags. 
  orf.nodes <- .Object@de.genes %>%
    dplyr::select(alt.name = name, name = ORF, logFC) %>%
    mutate(group = case_when(
      logFC < 0 ~ 'down',     #When logFC < 0 its down-regulated
      TRUE ~ 'up'             #In all other cases it's up (note that 0 is not possbile as these are DE genes only)
    )) %>%
    mutate(size = rep(1,n())) %>% 
    dplyr::select(-logFC) %>%
    unique()
  
  # - Defining the pathway nodes, and adjusting the size 
  # - according to the number of connected nodes
  pathway.nodes <- .Object@de.genes %>%
    group_by(pathway.name) %>%
    summarise(size = n()) %>%
    mutate(group = rep('pathway',n()),
           alt.name = pathway.name) %>% # just use same name as alternative name
    dplyr::select(name = pathway.name, alt.name, group, size)
  
  # - Combining the ORF and pathway nodes into a single data frame
  nodes <- rbind(orf.nodes, pathway.nodes)
  cat(paste0(nrow(nodes), ' nodes were created'))
  
  .Object@nodes <- nodes
  .Object
})

setMethod("connect.nodes", signature("KEGG.network"), function(.Object) {
  # - The edges in the network are defined by mentioning the index of the target 
  # - and source node in the `nodes` data frame, rather then specifing the node 
  # - names itself. In defining the edges (links) we need to know the indeces of 
  # - each node and therefore an linking table would make this much easier. 
  nodes.indexed <- .Object@nodes %>%
    mutate(index = seq(0, n()-1))
  
  # - building a link set; defining which nodes are connected.
  # - Important here is that the nodes are coupled based on 
  # - their position in the nodes data frame (i.e. the index)
  links <- data.frame(
    source <- .Object@de.genes %>%
      left_join(nodes.indexed, by = c('ORF'='name')) %>%
      .$index,
    target <- .Object@de.genes %>%
      left_join(nodes.indexed, by = c('pathway.name'='name')) %>%
      .$index
  )
  links$value <- rep(1, nrow(links))
  colnames(links) <- c('source','target','value')
  
  .Object@links <- links
  .Object
}
)


setMethod("plot", signature("KEGG.network"), function(.Object, name, direct.visible = FALSE) {
  opacity <- if (direct.visible) 1 else 0
  network <- forceNetwork(Links = .Object@links, Nodes = .Object@nodes, Source = "source",
                          Target = "target", Value = "value", NodeID = "name",
                          linkWidth = 2,  Nodesize = 'size', Group = "group", opacity = 0.9 , zoom = TRUE,
                          charge = -30, opacityNoHover = opacity, fontSize = 10, linkDistance = 40,
                          linkColour = "#000",  colourScale = 'd3.scaleOrdinal().domain(["up", "down", "pathway"]).range(["red", "green", "#8B008B"])')
  return (network)
}
)

