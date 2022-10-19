# Everything up to DESeq is the same

# Load libraries
library(tidyverse)
library(tximport)
library(DESeq2)
library(clusterProfiler)
library(WGCNA)
library(pheatmap)

setwd("~/Transcriptomic/julia/") # change to the directory with kallisto results

# Create a vector with sample names
sample_names <- c(paste0("Contr", 1:3), paste0("Pdrm", c(1:3)))
# Create a list of kallisto directories, make sure it corresponds to sample names in order and length
kallisto_dirs <- list.dirs(".")[-1]

# Create a phenotypic table
samples <- data.frame(sample = sample_names,
                      condition = c("Contr", "Contr","Contr","Pdrm","Pdrm","Pdrm"),
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
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "Contr")

# Calculate DE
dds <- DESeq(ddsTxi)
dds <- DESeq(dds)

# another way of doing the same
# dds <- ddsTxi %>% estimateSizeFactors %>% estimateDispersions %>% nbinomWaldTest

# Check normalization factors
normalizationFactors(dds)

resultsNames(dds)  # check contrasts names
# get results
res <- results(dds, name = "condition_Contr_vs_Pdrm")

# same way of getting results
# results(dds, contrast = c("time", "T6", "T3")) yt рабочая команда тк нет тайм

# Check model matrix
model.matrix(~+condition, samples) %>% View()
model.matrix(~0+condition, samples) %>% View()

# summary of results
summary(res)

# Gene counts for a specific gene
plotCounts(dds, gene="ENSMUST00000197754.1", intgroup="condition")

# Plot volcano
res %>%
  as.data.frame %>% 
  ggplot(aes(log2FoldChange, -log10(padj), color = padj < 0.05))+
  geom_point()+
  scale_color_manual(values=c("black", "red"))+
  xlim(c(-2.5, 2.5))


## Concordance tests

# Plot a heatmap of 100 most expressed genes
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("condition")])
vsd <- varianceStabilizingTransformation(dds)
pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=T,
         cluster_cols=T, annotation_col=df)

# Distance between samples heatmap
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$condition, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         color = colors)


# variance-stabilizing transformation of data (important for clustering or ML)


plotPCA(vsd, intgroup = c("condition"))

### Enrichment ###

BiocManager::install("org.Mm.eg.db")

# Filter only significant results
sign_results <- res %>%
  as.data.frame %>%
  rownames_to_column("gene_name") %>%
  mutate(ens_id = str_replace(gene_name, "\\.[1-9]+", "")) %>%
  filter(padj < .05)

# check up and down regulated genes separately
sign_up <- sign_results %>% filter(log2FoldChange > 0)
sign_dw <- sign_results %>% filter(log2FoldChange < 0)


# KEGG enrichment
KEGG_enrich <- enrichKEGG(sign_dw$ens_id, organism = "mus")

# vizualize results
barplot(KEGG_enrich)
dotplot(KEGG_enrich)
cnetplot(KEGG_enrich)

# let's see how DE genes map to the pathway
browseKEGG(KEGG_enrich, "mus00500")


# calculate GO enrichment

GO_enrich <- enrichGO(sign_dw$ens_id, "org.Mm.eg.db", keyType = "ENSEMBL", ont = "MF")
dotplot(GO_enrich, showCategory = 20)

GO_enrich@result %>% View  # table view of enrichment results

GO_enrich <- enrichplot::pairwise_termsim(GO_enrich)
emapplot(GO_enrich)

# How enrichment maps to GO DAG
goplot(GO_enrich)






















