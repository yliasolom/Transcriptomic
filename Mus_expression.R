# Everything up to DESeq is the same
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
# Load libraries
library(tidyverse)
library(tximport)
library("DESeq2")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rhdf5")
library(rhdf5)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

library("clusterProfiler")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("WGCNA")

library(WGCNA)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pheatmap")

library(pheatmap)

setwd("/home/aglab/projects/transcriptomics/ptj/GSE121884") # change to the directory with kallisto results

# Create a vector with sample names
sample_names <- c(paste0("WT", 1:3), paste0("KO", c(1:3)))

# Create a list of kallisto directories, make sure it corresponds to sample names in order and length
kallisto_dirs <- list.dirs(".")[-1]

# Create a phenotypic table
samples <- data.frame(sample = sample_names,
                      condition = c("WT", "WT","WT","KO","KO","KO"),
                      path = kallisto_dirs)

files <- file.path(samples$path, "abundance.h5")  # path to abundance files in kallisto dirs
txi <- tximport(files, type = 'kallisto', txOut = T)  # import transcripts

# Create DESeq object
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)  # change design if needed


# filter low-expressed genes
dds <- ddsTxi[rowSums(counts(ddsTxi)) >= 100,]

# relevel condition (time is fine cause T3 < T6)
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "WT")

# Calculate DE
dds <- DESeq(ddsTxi)
dds <- DESeq(dds)

# another way of doing the same
dds <- ddsTxi %>% estimateSizeFactors %>% estimateDispersions %>% nbinomWaldTest

# Check normalization factors
normalizationFactors(dds)

resultsNames(dds)  # check contrasts names
# get results
res <- results(dds, name = "condition_KO_vs_WT")

# same way of getting results
#results(dds, contrast = c("time", "T6", "T3"))

# Check model matrix
model.matrix(~+condition, samples) %>% View()
model.matrix(~0+condition, samples) %>% View()

# summary of results
summary(res)

# Gene counts for a specific gene
plotCounts(dds, gene="ENSMUST00000178537.1", intgroup="condition")

# Plot volcano
res %>%
  as.data.frame %>% 
  ggplot(aes(log2FoldChange, -log10(padj), color = padj < 0.05))+
  geom_point()+
  scale_color_manual(values=c("black", "red"))+
  xlim(c(-2.5, 2.5))

# variance-stabilizing transformation of data (important for clustering or ML)
vsd <- varianceStabilizingTransformation(dds)

plotPCA(vsd, intgroup = c("condition"))

## Concordance tests

# Plot a heatmap of 100 most expressed genes
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("condition")])
pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=T,
         cluster_cols=T, annotation_col=df)

# Distance between samples heatmap
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$time, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$condition, vsd$time, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         color = colors)




### Enrichment ###

BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

# Filter only significant results
sign_results <- res %>%
  as.data.frame %>%
  rownames_to_column("gene_name") %>%
  mutate(ens_id = str_replace(gene_name, "\\.[0-9]+", "")) %>%
  filter(padj < .05)

 


# check up and down regulated genes separately
sign_up <- sign_results %>% filter(log2FoldChange > 0)
sign_dw <- sign_results %>% filter(log2FoldChange < 0)
# не раюотает и не надо
kegg_res <- bitr(sign_results$ens_id, "ENSEMBLTRANS", "ENTREZID", OrgDb = "org.Mm.eg.db")
keggup_res <- bitr(sign_up$ens_id, "ENSEMBLTRANS", "ENTREZID", OrgDb = "org.Mm.eg.db")
# KEGG enrichment
KEGG_enrich <- enrichKEGG(kegg_res$ENTREZID, organism = "mmu")
KEGG_enrichup <- enrichKEGG(keggup_res$ENTREZID, organism = "mmu")
# vizualize results
barplot(KEGG_enrichup)
dotplot(KEGG_enrich)
cnetplot(KEGG_enrich)


####################################
# let's see how DE genes map to the pathway
browseKEGG(KEGG_enrich, "mmu00500")

BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

# calculate GO enrichment
GO_enrich <- enrichGO(sign_dw$ens_id, "org.Mm.eg.db", keyType = "ENSEMBLTRANS", ont = "MF")
dotplot(GO_enrich, showCategory = 20)

GO_enrich@result %>% View  # table view of enrichment results

GO_enrich <- enrichplot::pairwise_termsim(GO_enrich)
install.packages("ggnewscale")
emapplot(GO_enrich)

install.packages("ggnewscale")

# How enrichment maps to GO DAG
goplot(GO_enrich)

# Semantic similarity
scGO <- GOSemSim::godata("org.Mm.eg.db", ont = "MF")
# сравнение делать не надо 
# Compare two MF terms 
GOSemSim::goSim("GO:0005215", "GO:0022857", semData=scGO, measure="Jiang") погуглить цифирки хроматин билдинг и хелекасе фстивити 

# Compare two groups of terms
go1 = paste0("GO:", c("0020037", "0046906"))
go2 = paste0("GO:", c("0015002", "0016675"))


GOSemSim::mgoSim(go1, go2, semData=scGO, measure="Jiang", combine = "BMA")

# write this table to a file
sign_up %>% write_tsv("sign_up.tsv")
##########################################


# FGSEA
BiocManager::install("rWikiPathways")
rWikiPathways::downloadPathwayArchive(organism="Mus musculus", format="gmt")
# load pathways
pathways <- fgsea::gmtPathways("wikipathways-20221010-gmt-Mus_musculus.gmt")

# create ranks and names
ranks_for_gsea <- res %>% 
  as.data.frame() %>%   
  na.omit() %>% 
  arrange(desc(stat)) %>% 
  pull(stat)

names(ranks_for_gsea) <- res %>% 
  as.data.frame() %>% 
  na.omit() %>%
  rownames_to_column("gene_name") %>%
  mutate(ens_id = str_replace(gene_name, "\\.[0-9]+", "")) %>% 
  pull(ens_id) %>% 
  bitr("ENSEMBLTRANS", "ENTREZID", "org.Mm.eg.db") %>% 
  pull(ENTREZID)  ## change names to symbols

# run fgsea
fgsea_results <- fgsea::fgsea(pathways, ranks_for_gsea)

# inspect and plot results
head(fgsea_results[order(pval), ])

# GSEA plot for one pathway (here chosen arbitrarily)
fgsea::plotEnrichment(pathways[["mRNA processing%WikiPathways_20221010%WP310%Mus musculus"]], ranks_for_gsea) + labs(title="mRNAprocess")

# plot top 10 positively and negatively enriched gene lists
topPathwaysUp <- fgsea_results[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgsea_results[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
fgsea::plotGseaTable(pathways[topPathways], ranks_for_gsea, fgsea_results, 
                     gseaParam=0.5)


# clusterProfiler does fgsea too!
# it needs different names though

names(ranks_for_gsea) <- res %>% 
  as.data.frame() %>% 
  na.omit() %>%
  arrange(desc(stat)) %>% 
  rownames_to_column("gene_name") %>%
  mutate(ens_id = str_replace(gene_name, "\\.[0-9]+", "") %>% str_replace("transcript:", "")) %>%
  pull(ens_id) %>% 
  bitr("ENSEMBLTRANS", "ENTREZID", "org.Mm.eg.db") %>% 
  pull(ENTREZID)  ## change names to symbols

GO_gsea <- gseGO(ranks_for_gsea, ont = "ALL", org.Mm.eg.db, eps = 0) 

# install.packages("ggridges")
# save picture to a file. Plot has to be flanked by png and dev.off commands
png("ridge.png", 1600, 900)
install.packages("ggridges")
ridgeplot(GO_gsea)

dev.off()

# plot GSEA results for the first GO category
gseaplot(GO_gsea, 1, title = GO_gsea@result$Description[[1]])


### WGCNA ###
#не надо делать
# Expression data (has to be transformed)
datExpr <- t(assay(vsd))
rownames(datExpr) <- samples$sample

# read metabolome data
datTraits <- read_tsv("table.tsv") %>% 
  dplyr::select(-`8.mzXML`) %>% 
  column_to_rownames('name')

datTraits <- cbind(datTraits, datTraits)  # double the table because we decided to use SE reads

datTraits <- t(datTraits)

traitColors <- numbers2colors(datTraits, signed = FALSE)


# Clustering samples
sampletree <- hclust(dist(datExpr), method = "average")
plot(sampletree, cex = 0.5)

# Plot how samples and phenotypes correspond
plotDendroAndColors(sampletree, traitColors,
                    groupLabels = colnames(datTraits), cex.colorLabels = 0.5,
                    main = "Sample dendrogram and trait heatmap")

# Pick powers for coexpression network construction
powers <- c(c(1:10), seq(from = 15, to=50, by=5))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)



# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# construct the network
net <- blockwiseModules(datExpr, power = 8,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = FALSE,
                        saveTOMFileBase = "yeastTOM",
                        verbose = 3)



module_colors <- labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)

# Extract module eigengenes
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 12)


# Plot module behavior across the dataset
MEs %>%
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  ggplot(aes(sample, MEblue))+
  geom_point()+
  theme_bw()



# Plot module-phenotype correlation
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# analyze some modules
blue_module <- colnames(datExpr)[module_colors=="blue"] %>% str_replace("_mRNA", "")
red_module <- colnames(datExpr)[module_colors=="red"] %>% str_replace("_mRNA", "")

# Are modules enriched in something?
KEGG_enrich <- enrichKEGG(blue_module, organism = "sce")


dotplot(KEGG_enrich)

GO_enrich <- enrichGO(blue_module, "org.Sc.sgd.db", keyType = "ENSEMBL", ont = "ALL")
barplot(GO_enrich)


# we can also visualize network data in Cytoscape. Created files can be imported there
datexpr_brown <- datExpr[, module_colors=="green"]
TOM_brown = TOMsimilarityFromExpr(datexpr_brown, power = 8, networkType = "signed", TOMType="signed")
probes = colnames(datexpr_brown)
dimnames(TOM_brown) = list(probes, probes)
exportNetworkToCytoscape(TOM_brown, "edges.txt", "nodes.txt", nodeNames = probes, threshold = 0.1)




