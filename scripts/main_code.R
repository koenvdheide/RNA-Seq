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
library(tidyr)
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
source('scripts/RNA_seq_analyzer.R')
source('scripts/KEGG_annotater.R')
source('scripts/KEGG_network.R')

##################################################################################################
#------------------------------------------------------------------------------------------------
# PART 1: Loading data + KEGG analysis
#------------------------------------------------------------------------------------------------
##################################################################################################

##################################################################################################
# Loading data
##################################################################################################
SPECIES = 'lpl'

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
# - Swith to output folder
setwd('/home/rick/Desktop/HAN/jaar 3/HAN/RNA_seq_clean/output')
write.table(wcfs1.de, file = 'wcfs1_de_genes.txt',sep = '\t', row.names = F)

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

# - Running the KEGG annotation code
wcfs1.kegg <- new('KEGG.annotater',wcfs1.result, wcfs1.de, 'lpl')
wcfs1.kegg <- gather.kegg.annotation(wcfs1.kegg)
wcfs1.kegg <- couple.kegg.annotation(wcfs1.kegg)
wcfs1.kegg <- filter.kegg.pathways(wcfs1.kegg, priority.interest)

# - Looking at the pathways we found and their ids
kegg.table <- get.kegg.table(wcfs1.kegg)
print(kegg.table)

# - We can download all pathways at once
map.kegg.pathway(wcfs1.kegg, kegg.table$pathway.id, kegg.table$pathway.name)

# - Or specifcy specific pathways (such as the one of our interest)
wanted <- kegg.table %>% filter(pathway.name %in% priority.interest)
map.kegg.pathway(wcfs1.kegg, c('lpl03010','lpl03018'), c('Ribosome','RNA degradation'))

# - As pathview does not automatically remove the pathway PNG and XML files
# - We wrote a little function to do so
remove.junk(wcfs1.kegg, kegg.table$pathway.id )


##################################################################################################
# WCFS1 analysis: Plotting the fold changes for the pathways
##################################################################################################
plot.fold.changes(wcfs1.kegg, type = 1)
plot.fold.changes(wcfs1.kegg, type = 2)

##################################################################################################
# WCFS1 analysis: saving DE + KEGG annotation to file
##################################################################################################
write.table(wcfs1.kegg@de.kegg.genes, file = 'wcfs1_de_genes_kegg.txt',sep = '\t', row.names = F)

##################################################################################################
# WCFS1 analysis: KEGG GSEA network using own code 
##################################################################################################
de.genes.annotated <- wcfs1.kegg@de.kegg.genes.filtered %>% left_join(anno)
wcfs1.network <- new('KEGG.network',de.genes.annotated)
wcfs1.network <- build.nodes(wcfs1.network)
wcfs1.network <- connect.nodes(wcfs1.network)

# Show network in which user needs to move over the nodes to see the name
network <- draw.network(wcfs1.network, 'name')
network

# Showing the genes names in the network, so hovering is not needed
network.named <- draw.network(wcfs1.network, 'alt.name', direct.visible =  TRUE)
network.named

##################################################################################################
# WCFS1 analysis: saving the network
##################################################################################################
saveNetwork(network.named, file = 'wcfs1_network.html', selfcontained = TRUE)

#################################################################################################
#------------------------------------------------------------------------------------------------
# PART 2: KEGG analysis using ClusterProfiler package
#------------------------------------------------------------------------------------------------
#################################################################################################
# - This package did not really fit our needs, hence we wrote our own code (see above). We 
# - decided to still include this code at the bottom as some of its plots were useful and it was
# - still nice to have a comparison. 

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
                     minGSSize = 1) #minimal size of each geneSet for analyzing

# - Plotting the fold changes of the genes and their connection to the pathway
# - in which these were found. 
cnet.plot <- cnetplot(kegg.gsea, foldChange = wcfs1.gene.list, showCategory = priority.interest,
              node_label = FALSE)

# Hacky way to remove gene labels and increase size of pathway labels
cnet.plot$data$name <- gsub("lp.*"," ",cnet.plot$data$name)
cnet.plot + geom_node_text(aes_(label=~name), size = 5)

# Additional KEGG gsea plots
emapplot(kegg.gsea, showCategory = priority.interest)
heatplot(kegg.gsea, foldChange=wcfs1.gene.list)
ridgeplot(kegg.gsea)

# Saving the rigdeplot for the poster
png("wcf1_kegg_gsea_ridge_plot.png", width = 30, height = 15, units = 'cm', res = 400)
ridgeplot(kegg.gsea)
dev.off()

##################################################################################################
# WCFS1 analysis: Why is the pentose phosphate pathway not enriched by ClusterProfiler?
##################################################################################################
# Selecting all genes and look their pathways and color the according to whether these 
# were significant (DE)
wcfs1.sign.table <- wcfs1.kegg@all.kegg.genes.filtered %>%
  mutate(significant = 
           case_when(
             adjpvalue < 0.05 ~ TRUE,
             TRUE ~ FALSE # adjpvalue >= 0.05 to FALSE
           ))

# Plotting the number of signifcant genes and non significant genes for 
# each of the pathways of interest
wcfs1.sign.table %>%
  group_by(pathway.name, significant) %>%
  summarise(count = n()) %>% # Counting the number of signifcant and non-significant genes per pathway
  ggplot(aes( x = pathway.name, y = count, fill = significant)) +
  geom_bar(stat = 'identity')


# Plotting the fold change spread over the pahtways and color
# according to whether the genes was signifcant or not
wcfs1.sign.table %>%
  ggplot(aes(x = pathway.name, y = logFC, col = significant )) +
  geom_point(size = 3)
















