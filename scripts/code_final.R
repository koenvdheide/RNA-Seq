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
source('scripts/RNA_seq_analyzer_class.R')
source('scripts/kegg_annotater_class.R')

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
wcfs1 <- new('RNA.seq.analyzer',count.data[,1:4], anno)
wcfs1 <- filter.low.counts(wcfs1)
wcfs1 <- run.all(wcfs1)
wcfs1 <- determine.significance(wcfs1, 0.05)

nc8 <- new('RNA.seq.analyzer',count.data[,5:8], anno)
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

# - This will return the original edgeR DGEList object, so it can be passed to all edgeR functions easily
wcfs1.edgerobj <- wcfs1@edger.data 

# - This will give the differentially expressed (DE) genes as identified by edgeR
wcfs1.de <- wcfs1@edger.de %>%   rownames_to_column('ORF')

# - This will give all genes (both DE and not DE)
wcfs1.result <- wcfs1@edger.result %>%   rownames_to_column('ORF')

##################################################################################################
# WCFS1 analysis: saving DE genes to a file
##################################################################################################
write.table(wcfs1.de, file = 'output/wcfs1_de_genes.txt',sep = '\t', row.names = F)

##################################################################################################
# WCFS1 analysis: MA plot
##################################################################################################
plotSmear(wcfs1.edgerobj, de.tags = wcfs1.de$ORF)

##################################################################################################
# WCFS1 analysis: Gather KEGG annotation and filtering it
##################################################################################################
# - As one gene can correspond to multiple pathways we will specifcy our pathways of interest
priority.interest <- c('Pyruvate metabolism','Glycolysis / Gluconeogenesis',
                       'Fatty acid biosynthesis','Galactose metabolism',
                       'Pentose phosphate pathway','Pyrimidine metabolism')

wcfs1.kegg <- new('KEGG.annotater',wcfs1.result, wcfs1.de, 'lpl')
wcfs1.kegg <- get.kegg.annotation(wcfs1.kegg)
wcfs1.kegg <- couple.kegg.annotation(wcfs1.kegg)
wcfs1.kegg <- filter.kegg.pathways(wcfs1.kegg, priority.interest)
map.kegg.pathway(wcfs1.kegg, c('lpl00010'))


##################################################################################################
# WCFS1 analysis: Plotting the fold changes for the pathways
##################################################################################################
plot.fold.changes(wcfs1.kegg, type = 1)
plot.fold.changes(wcfs1.kegg, type = 2)

##################################################################################################
# WCFS1 analysis: saving DE + KEGG annotation to file
##################################################################################################
write.table(wcfs1.kegg@de.kegg.genes, file = 'output/wcfs1_de_genes_kegg.txt',sep = '\t', row.names = F)

##################################################################################################
# WCFS1 analysis: KEGG GSEA network using own code 
##################################################################################################
de.genes.annotated <- wcfs1.kegg@de.kegg.genes.filtered %>% left_join(anno)
wcfs1.network <- new('KEGG.network',de.genes.annotated)
wcfs1.network <- build.nodes(wcfs1.network)
wcfs1.network <- connect.nodes(wcfs1.network)

# Show network in which user needs to move over the nodes to see the name
network <- plot(wcfs1.network, 'name')
network

# Showing the genes names in the network, so hovering is not needed
network.named <- plot(wcfs1.network, 'alt.name', TRUE)
network.named

##################################################################################################
# WCFS1 analysis: Saving the network
##################################################################################################
original.wd <- getwd()
setwd(paste0(original.wd,'/output'))
saveNetwork(network.named, file = 'network_ribose_glucose.html', selfcontained = TRUE)

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
dev.off(

















