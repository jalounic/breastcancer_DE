dev.off()
rm(list = ls())

#library(BiocManager)
library(Rsubread)
library(DESeq2)
library(stringr)
library(pheatmap)
library(clusterProfiler)
library(tidyverse)


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
drop <- c("Chr", "Start", "End", "Strand", "Length")
fC_trimmed <- fC_out[!names(fC_out) %in% drop]

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
# put into function, do for all 3 comparisons?

NormvsTNBC <- results(dds, contrast=c("condition","Normal", "TNBC"))
TNBCvsNorm <- results(dds, contrast=c("condition", "TNBC", "Normal"))

#
plotCounts(dds, "ENSG00000285701", intgroup = "condition", normalized = TRUE)

# How many genes are differentially expressed (DE) in the pairwise comparison you selected (e.g. padj < 0.05)

# NAs:
ldata <- as.data.frame(NormvsTNBC@listData)
ldata <- cbind(NormvsTNBC@rownames, ldata)

DEgenes <- ldata[which(ldata$padj < 0.05),]

colnames(DEgenes)[1] <- "geneid"
colnames(DEgenes)


DEgenes <- DEgenes %>% 
  mutate(upregulated = log2FoldChange < 0)



results(dds, contrast=c("condition","Normal", "NonTNBC"))
results(dds, contrast=c("condition","NonTNBC", "TNBC"))










