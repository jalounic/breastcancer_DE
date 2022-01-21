dev.off()
rm(list = ls())

#library(BiocManager)
library(Rsubread)
library(DESeq2)
library(stringr)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(enrichplot)

# packageVersion("Rsubread")

# import featurecounts table:
fC_out <- read.delim("~/Desktop/fC_out", comment.char="#", row.names = "Geneid")

#clean up column names
colnames(fC_out)

splitted <- str_split(names(fC_out)[str_detect(names(fC_out), "X.")], '\\.')

samples <- c()
for (line in splitted){
  samples <- c(samples, line[9])
}
names(fC_out)[str_detect(names(fC_out), "X.")] <- samples
colnames(fC_out) # much better!

# drop columns:
fC_trimmed <- fC_out[!names(fC_out) %in% c("Chr", "Start", "End", "Strand", "Length")]

conditions <- c(rep("NonTNBC", 3), rep("Normal", 3), rep("TNBC", 3))

# make DESed df!
dseqobj <- DESeqDataSetFromMatrix(
  countData = fC_trimmed,
  design = ~ condition,
  colData = data.frame(
    condition = conditions))

dds <- DESeq(dseqobj)
vsd <- vst(dds, blind=TRUE) # remove dependence of variance to mean

plotPCA(vsd)

## heatmaps!!! (include in report or no?)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:30]
df <- as.data.frame(colData(dds)["condition"])
rownames(df) <- colnames(dds)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# now let's get results!

#look at potentially interesting genes:
# ENSG00000062038 (CDH3)
plotCounts(dds, "ENSG00000062038", intgroup = "condition", normalized = TRUE) # upregulated in TNBC, down in NonTNBC!

# ENSG00000128422 (KRT17)
plotCounts(dds, "ENSG00000128422", intgroup = "condition", normalized = TRUE) # upregulated in TNBC, down in NonTNBC!

# How many genes are differentially expressed (DE) in the pairwise comparison you selected (e.g. padj < 0.05)

getDEGenes <- function(deseq2results){
  ldata <- as.data.frame(deseq2results@listData) #60'664 total observations
  ldata <- cbind(deseq2results@rownames, ldata)
  
  DEgenes <- ldata[which(ldata$padj < 0.05),]
  
  colnames(DEgenes)[1] <- "geneid"
  colnames(DEgenes)
  
  DEgenes <- DEgenes %>% 
    mutate(upregulated = log2FoldChange < 0)
  DEgenes <- DEgenes %>% 
    mutate(downregulated = log2FoldChange > 0)
  return(DEgenes)
}

# look for specific genes: 
# DEgenes$log2FoldChange[DEgenes$geneid == "ENSG00000128422"]

NormvsTNBC <- getDEGenes(results(dds, contrast=c("condition","Normal", "TNBC")))

sum(NormvsTNBC$upregulated) # 13'698
sum(NormvsTNBC$downregulated) # 4'973
nrow(NormvsTNBC) # 18'671

#look whether genes of interest are in there
"ENSG00000128422" %in% NormvsTNBC$geneid

NormvsNonTNBC <-getDEGenes(results(dds, contrast=c("condition","Normal", "NonTNBC")))

sum(NormvsNonTNBC$upregulated) # 10'825
sum(NormvsNonTNBC$downregulated) # 4'681
nrow(NormvsNonTNBC) # 15'506

#look whether genes of interest are in there
"ENSG00000128422" %in% NormvsNonTNBC$geneid


NonTNBCvsTNBC <- getDEGenes(results(dds, contrast=c("condition","NonTNBC", "TNBC")))

sum(NonTNBCvsTNBC$upregulated) # 1172 (TNBC > NonTNBC)
sum(NonTNBCvsTNBC$downregulated) # 1210 (NonTNBC > TNBC)
nrow(NonTNBCvsTNBC) # 2382 total differentially expressed genes

#look whether genes of interest are in there
"ENSG00000128422" %in% NonTNBCvsTNBC$geneid

## now let's get those GO terms!
DE_IDs <- union(NormvsTNBC$geneid, NormvsNonTNBC$geneid)
DE_IDs <- union(DE_IDs, NonTNBCvsTNBC$geneid)
length(DE_IDs) # 22'157 total DE genes

universe <- rownames(fC_out) # 60'664

ego <- enrichGO(gene          = DE_IDs,
                universe      = universe,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE,
                pool = TRUE)
head(ego)

?clusterProfiler::enrichGO


barplot(ego, showCategory=10) 
dotplot(ego, showCategory=10)

## convert gene ID to Symbol
p1 <- cnetplot(ego, showCategory = 5, node_label="category", 
               cex_label_category = 1.2) 
p2 <- cnetplot(ego, node_label="gene", 
               cex_label_gene = 0.8) 
p3 <- cnetplot(ego, node_label="all") 
p4 <- cnetplot(ego, node_label="none", 
               color_category='firebrick', 
               color_gene='steelblue') 
p1
p2
p3
p4
