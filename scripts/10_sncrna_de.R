#!/usr/bin/env Rscript
# GSE332551 — small-RNA DESeq2 (AD vs CN), classified, for miRNA-mRNA integration.
# Writes results/sncrnaseq_DE_all.csv (all tested features + class + direction).
suppressMessages(library(DESeq2))

counts <- read.csv("data/GSE332551_counts_sncRNAseq.csv", row.names = 1, check.names = FALSE)
libsize <- colSums(counts)
cat("Small-RNA library sizes:\n"); print(libsize)

cond <- factor(ifelse(grepl("AD$", colnames(counts)), "AD", "CN"), levels = c("CN", "AD"))
coldata <- data.frame(condition = cond, row.names = colnames(counts))

keep <- rowMeans(counts) >= 5
dds <- DESeqDataSetFromMatrix(round(as.matrix(counts[keep, ])), coldata, ~condition)
dds <- DESeq(dds, quiet = TRUE)
cat("\nDESeq2 size factors:\n"); print(round(sizeFactors(dds), 3))

res <- as.data.frame(results(dds, contrast = c("condition", "AD", "CN")))
res$sncRNA <- rownames(res)

classify <- function(x) {
  ifelse(grepl("^hsa-(miR|let)", x), "miRNA",
  ifelse(grepl("^hsa-piR", x), "piRNA",
  ifelse(grepl("^tRF", x, ignore.case = TRUE), "tRF",
  ifelse(grepl("SNOR|^U[0-9]|snoR", x, ignore.case = TRUE), "snoRNA", "other"))))
}
res$class <- classify(res$sncRNA)
res$direction <- ifelse(res$log2FoldChange >= 0, "up_in_AD", "down_in_AD")

res <- res[order(res$pvalue), c("sncRNA", "class", "baseMean", "log2FoldChange",
                                "lfcSE", "pvalue", "padj", "direction")]
write.csv(res, "results/sncrnaseq_DE_all.csv", row.names = FALSE)

de <- subset(res, !is.na(pvalue) & pvalue < 0.05 & abs(log2FoldChange) > 1)
cat(sprintf("\nTested %d features; %d DE (p<0.05, |log2FC|>1)\n", nrow(res), nrow(de)))
cat("DE miRNAs by direction:\n")
print(table(de$class, de$direction))
