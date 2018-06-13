##################################################################################################
# Authors: Rick Beeloo, Koen van der Heide, and Thomas Reinders
# Summary: This pipeline processes RNA seq counts:
# - Filtering low counts
# - Normalizing the counts
# - Showing a package comparison (DESeq2, EdgeR, and Limma)
# - Coupling the genes to KEGG pathways 
# - Performing KEGG analysis (such as GSEA, network, FC spread)
##################################################################################################

# Data formatting libraries
library(tibble)
library(dplyr)

# Plotting libraries
library(ggplot2)
library(ggraph)
library(VennDiagram)
library(enrichplot) # <- install from github
library(networkD3)
library(pathview)

# RNA specific libraries
library(DESeq2)
library(edgeR)
library(limma)

# Annotation libraries
library(clusterProfiler)

# Setting working dir
setwd('/home/rick/Desktop/HAN/jaar 3/HAN/RNA_seq_clean')

# Loading the other scripts we need
source('scripts/venn.extract.R')
source('scripts/class.R')

##################################################################################################
# GLOBAL VARIABLES
##################################################################################################
SPECIES = 'lpl'

##################################################################################################
# Loading data
##################################################################################################
# - RNA seq counts
count.data <- as.matrix(read.table('data/RNA-Seq-counts.txt', sep = '\t', header = T, skip = 1, row.names = 1))

# - Outdated annotation 
anno <- read.table('data/WCFS1_anno.txt', sep = '\t', header = T, fill = T, quote = "")

# - Add Entrez ids 
anno <- anno %>% left_join(bitr_kegg(anno$ORF, fromType='kegg', toType='ncbi-geneid', organism=SPECIES), by = c('ORF'='kegg'))


##################################################################################################
# Running limma, edgeR and DESeq2
##################################################################################################
wcfs1 <- new('SeqDat',count.data[,1:4], anno)
wcfs1 <- filter.low.counts(wcfs1)
wcfs1 <- run.all(wcfs1)
wcfs1 <- determine.significance(wcfs1, 0.05)

nc8 <- new('SeqDat',count.data[,5:8], anno)
nc8 <- filter.low.counts(nc8)
nc8 <- run.all(nc8)
nc8 <- determine.significance(nc8, 0.05)

##################################################################################################
# Differences between RNA-seq packages
##################################################################################################
wcfs1 <- draw.venn.diagram(wcfs1)
nc8 <- draw.venn.diagram(nc8)

##################################################################################################
# WCFS1 analysis: getting edgeR data
##################################################################################################
# - We choose to use the edgeR data for further analysis
# - If the core set is prefered use wcfs1@venn.data$limmaedgerdeseq2 

# - This will return the original edgeR DGEList object, so it can be passed to all edgeR functions
# - easily
wcfs1.edgerobj <- wcfs1@edger.data 

# - This will give the differentially expressed (DE) genes as identified by edgeR
wcfs1.de <- wcfs1@edger.de %>%
  rownames_to_column('ORF')

# - This will give all genes (both DE and not DE)
wcfs1.result <- wcfs1@edger.result %>%
  rownames_to_column('ORF')

##################################################################################################
# WCFS1 analysis: saving DE genes to a file
##################################################################################################
write.table(wcfs1.de, file = 'output/wcfs1_de_genes.txt',sep = '\t', row.names = F)

##################################################################################################
# WCFS1 analysis: MA plot
##################################################################################################
plotSmear(wcfs1.edgerobj, de.tags = wcfs1.de$ORF)

# For the poster we need to increase the font size and image size
png("output/wcfs1_ma_plot.png", width = 30, height = 15, units = 'cm', res = 400)
plotSmear(wcfs1.edgerobj, de.tags = wcfs1.de$ORF, cex.axis = 2, cex.lab = 1.5, cex = 0.8)
dev.off()

##################################################################################################
# WCFS1 analysis: GSEA KEGG using ClusterProfiler
##################################################################################################
# - Format data for Enrichment analsyis
# - Making a ranked list of genes based on their fold changes
# - Note that we use all genes here not the DE only
wcfs1.gene.list <- wcfs1.result$logFC
names(wcfs1.gene.list) <- wcfs1.result$ORF
wcfs1.gene.list = sort(wcfs1.gene.list, decreasing = TRUE)

# - Conducting KEGG GSEA using the ClusterProfiler package
kegg.gsea <- gseKEGG(geneList  = wcfs1.gene.list,
                     organism = SPECIES,
                     pvalueCutoff = 0.05, 
                     minGSSize = 1)

# - As one gene can correspond to multiple pathways we will use the previously
# - conducted KEGG enrichment + knowledge about ribose pathways to specify a priority/selection
priority.interest <- c('Pyruvate metabolism','Glycolysis / Gluconeogenesis',
                       'Fatty acid biosynthesis','Galactose metabolism',
                       'Pentose phosphate pathway','Pyrimidine metabolism')

# - Plotting the fold changes of the genes and their connection to the pathway
# - in which these were found. 
p <- cnetplot(kegg.gsea, foldChange = wcfs1.gene.list, showCategory = priority.interest,
              node_label = FALSE)

# Hacky way to remove gene labels and increase size of pathway labels
p$data$name <- gsub("lp.*"," ",p$data$name)
p + geom_node_text(aes_(label=~name), size = 5)

# Additional KEGG gsea plots
emapplot(kegg.gsea, showCategory = priority.interest)
heatplot(kegg.gsea, foldChange=wcfs1.gene.list)
ridgeplot(kegg.gsea)

# Saving the rigdeplot for the poster
png("output/wcf1_kegg_gsea_ridge_plot.png", width = 30, height = 15, units = 'cm', res = 400)
ridgeplot(kegg.gsea)
dev.off()

##################################################################################################
# WCFS1 analysis: Coupling genes to KEGG
##################################################################################################
# - Although the ClusterProfiler package produces nice plots (see code above) it 
# - does not allow for any flexibiltity in the plot it creates. Therefore we decided
# - to retrieve KEGG data our selfs and conduct our analysis and plotting

# - Retrieving data from KEGG and make annotation table hereof. Note that
# - we could have used the KEGGREST package. However this package
# - returns a character vector insted of a data frame, which would require 
# - an extra processing step. 
pathway.id.name <- read.table(paste0("http://rest.kegg.jp/list/pathway/",SPECIES), quote="", sep="\t")
gene.pathway.id <- read.table(paste0("http://rest.kegg.jp/link/pathway/",SPECIES), quote="", sep="\t")
gene.info <- left_join(gene.pathway.id, pathway.id.name, by = c('V2'='V1'))

# - Some formatting (makes downstream analysis easier)
colnames(gene.info) <- c("ORF",'pathway.id','pathway.name')
gene.info$ORF <- gsub('lpl:','',gene.info$ORF)
gene.info$pathway.id <- gsub('path:','',gene.info$pathway.id)
gene.info$pathway.name <- gsub(' - Lactobacillus plantarum WCFS1',"", gene.info$pathway.name)

# - Coupling KEGG annotation to the DE genes
wcfs1.de.kegg <- wcfs1.de %>% left_join(gene.info)

# - Let's see how many genes are coupled to a pathway
couple.count <- wcfs1.de.kegg %>%
  group_by(ORF) %>%
  summarise(coupled = !any(is.na(pathway.name))) %>%
  group_by(coupled) %>%
  summarise(number = n())

cat(paste0('Out of the ', n_distinct(wcfs1.de.kegg$ORF),' genes ',
           couple.count[1,2],' could be coupled to a KEGG entry, whereas ', 
           couple.count[2,2],' could not.'))

# - We filterd for pathways of interest to enhance interpretation as well as
# - comprehensiveness
wcfs1.de.kegg.filtered <- wcfs1.de.kegg %>%
  filter(pathway.name  %in% priority.interest) %>%
  arrange(-logFC)

# - A quick look at the genes and pathways
ggplot(wcfs1.de.kegg.filtered, aes(x = gsub(' ','\n',pathway.name), y=logFC, col = pathway.name)) + 
  geom_point(size = 2.5, show.legend = F) +
  xlab('Pathway') +
  ylab('LogFC') +
  ggtitle('Global pathways changes') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, col = 'red', size = 0.9) 

# - Alternatively, we could plot these genes ordered by FC and color according to the pathway
ggplot(wcfs1.de.kegg.filtered, aes(x = seq(1,nrow(wcfs1.de.kegg.filtered)), y = logFC, col = pathway.name)) +
  geom_point(size = 2.5)

##################################################################################################
# WCFS1 analysis: saving DE + KEGG annotation to file
##################################################################################################
write.table(wcfs1.de.kegg, file = 'output/wcfs1_de_genes_kegg.txt',sep = '\t', row.names = F)

##################################################################################################
# WCFS1 analysis: KEGG GSEA network using own code 
##################################################################################################
# As stated previously ClusterProfiler produces a nice network when using the 
# cnetplot command, however does not allow for flexibility in changing font sizes,
# node sizes, colors, spreading etc. Hene, we wrote our own code to produce a similar
# network. Besides the great variety of parameters that can be changed, the used package
# (networkD3) uses D3 javascript which allows (although indirectly) the export to SVG
# which can be formatted and adapted easily in adobe Illustrator for a poster.

# - Building an ORF node set, we colored these nodes based on whether 
# - the gene was identified as up or down regulated. Further 
# - The annotation is added as it more intuitive to interpret gene
# - names rather than the locus tags. 
orf.nodes <- wcfs1.de.kegg.filtered%>%
  left_join(anno) %>%
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
pathway.nodes <- wcfs1.de.kegg.filtered %>%
  group_by(pathway.name) %>%
  summarise(size = n()) %>%
  mutate(group = rep('pathway',n()),
         alt.name = pathway.name) %>% # just use same name as alternative name
  dplyr::select(name = pathway.name, alt.name, group, size)

# - Combining the ORF and pathway nodes into a single data frame
nodes <- rbind(orf.nodes, pathway.nodes)

# - The edges in the network are defined by mentioning the index of the target 
# - and source node in the `nodes` data frame, rather then specifing the node 
# - names itself. In defining the edges (links) we need to know the indeces of 
# - each node and therefore an linking table would make this much easier. 
nodes.indexed <- nodes %>%
  mutate(index = seq(0, n()-1))

# - building a link set; defining which nodes are connected.
# - Important here is that the nodes are coupled based on 
# - their position in the nodes data frame (i.e. the index)
links <- data.frame(
  source <- wcfs1.de.kegg.filtered %>%
    left_join(nodes.indexed, by = c('ORF'='name')) %>%
    .$index,
  target <- 
    wcfs1.de.kegg.filtered %>%
    left_join(nodes.indexed, by = c('pathway.name'='name')) %>%
    .$index
)
links$value <- rep(1, nrow(links))
colnames(links) <- c('source','target','value')

# - Optional to increase node size by a common factor
nodes$size <- nodes$size * 5

# - Quick function to making switching between ORF codes and ORF names cleaner
# - name: "name" will give ORF id, "alt.name" will give real gene name (if possible)
# - direct.visible: if FALSE hovering over the node will show the name, if TRUE
# -                 the node names will be directly visible. 
build.network <- function(name, direct.visible = FALSE) {
  opacity <- if (direct.visible) 1 else 0
  network <- forceNetwork(Links = links, Nodes = nodes, Source = "source",
                          Target = "target", Value = "value", NodeID = "name",
                          linkWidth = 2,  Nodesize = 'size', Group = "group", opacity = 0.9 , zoom = TRUE,
                          charge = -30, opacityNoHover = opacity, fontSize = 10, linkDistance = 40,
                          linkColour = "#000",  colourScale = 'd3.scaleOrdinal().domain(["up", "down", "pathway"]).range(["red", "green", "#8B008B"])')
  return (network)
}

# - Plot network with ORF codes
orf.code.net <- build.network('name')
orf.code.net

# - Plot network with ORF names
orf.name.net <- build.network('alt.name', TRUE)
orf.name.net

# - Saving the network, note that selfcontained  = TRUE will add all the required data 
# - for drawing inside the html file. 
# - As html widget can't be saved that easily we first have to change working directory
original.wd <- getwd()
setwd(paste0(original.wd,'/output'))
saveNetwork(orf.name.net, file = 'network_ribose_glucose.html', selfcontained = TRUE)

##################################################################################################
# WCFS1 analysis: Zooming in on the kegg pathways
##################################################################################################
# - After looking at the global changes regarding KEGG pathways it's useful to 
# - zoom in on these. For this task we used the commonly used PathView package.

# - Format the wcfs1.de data such that it can be handled by PathView
kegg.format <- wcfs1@edger.de %>%
  dplyr::select(logFC)

# - Pathview uses two files to create the output image:
# - PNG image of the KEGG pathay
# - XML file of the KEGG pathway
# - Those aren't deleted automatically, therefore we 
# - wrote this function
remove.junk <- function(path.id) {
  png.base = paste0(path.id, '.png')
  xml.file = paste0(path.id, '.xml')
  if (file.exists(png.base)) file.remove(png.base)
  if (file.exists(xml.file)) file.remove(xml.file)
}

# - This function will save the network for a given path id and 
# - output suffix. 
save.pathway <- function(path.id, out.suffix) {
  pathview(gene.data = kegg.format, 
           gene.idtype = "KEGG" , 
           pathway.id = path.id , 
           species = SPECIES,
           out.suffix = out.suffix, 
           map.symbol = TRUE,        # show gene names when possible
           same.layer = FALSE)       # extra layer for the gene names
  remove.junk(path.id)
}

# - To download a specific pathway
save.pathway('lpl00010', 'glycolysis')

# - Downloading all of them:
# - We already coupled the WCFS1 genes to each pathway so we can simply use this
# - same list to retrieve all pathway plots
pathway.info <-gene.info %>% 
  dplyr::select(pathway.id, pathway.name) %>%
  mutate(pathway.name = str_remove_all(.$pathway.name,'/')) %>% #remove / as pathview can't handle these when saving
  unique()

# - Gathering the data for each of the pathways in pathway.info
mapply(save.pathway, pathway.info$pathway.id, pathway.info$pathway.name)

# Move out of the output folder and back to the main folder
setwd(original.wd)

##################################################################################################
# WCFS1 analysis: Why is the pentose phosphate pathway not enriched by ClusterProfiler
##################################################################################################
# Note that we did this previously, however this time we use ALL the edgeR results,
# not only the DE ones as all the genes were used as input for ClusterProfiler.
wcfs1.sign.talbe <- wcfs1.result %>%
  dplyr::select(ORF,logFC, adjpvalue) %>%
  left_join(gene.info) %>%
  filter(pathway.name %in% priority.interest) %>%
  mutate(significant = 
           case_when(
             adjpvalue < 0.05 ~ TRUE,
             TRUE ~ FALSE # adjpvalue >-0.05 to FALSE
           ))

# Plotting the number of signifcant genes and non significant genes
wcfs1.sign.talbe %>%
  group_by(pathway.name, significant) %>%
  summarise(count = n()) %>%
  ggplot(aes( x = pathway.name, y = count, fill = significant)) +
  geom_bar(stat = 'identity')


# Plotting the fold change spread over the pahtways and color
# according to whether the genes was signifcant or not
wcfs1.sign.talbe %>%
  ggplot(aes(x = pathway.name, y = logFC, col = significant )) +
  geom_point(size = 3)





































