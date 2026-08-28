#!/usr/bin/env Rscript
# Diagnose why reanalysis DEGs (827, 3.6:1 up-skew) differ from paper (195, balanced)
suppressMessages(library(DESeq2))

counts_raw <- read.csv("../data/GSE332551_counts_RNAseq.csv", row.names = 1)
gene_names <- sub("^[^:]+:", "", rownames(counts_raw))
counts_gene <- aggregate(counts_raw, by = list(Gene = gene_names), FUN = sum)
rownames(counts_gene) <- counts_gene$Gene; counts_gene$Gene <- NULL

coldata <- data.frame(
  condition = factor(ifelse(grepl("AD", colnames(counts_gene)), "AD", "CN"),
                     levels = c("CN", "AD")),
  row.names = colnames(counts_gene))

cat("=== Library sizes (total counts per sample) ===\n")
print(colSums(counts_gene))

dds <- DESeqDataSetFromMatrix(round(counts_gene), coldata, design = ~ condition)
dds <- dds[rowMeans(counts(dds)) >= 10, ]
dds <- DESeq(dds, quiet = TRUE)

cat("\n=== DESeq2 size factors ===\n"); print(sizeFactors(dds))

# Sample clustering to spot outliers
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
cat("\n=== Sample correlation (Spearman, on VST) ===\n")
print(round(cor(assay(vsd), method = "spearman"), 3))
pc <- prcomp(t(assay(vsd)))
cat("\n=== PCA (PC1/PC2) ===\n")
print(round(pc$x[, 1:2], 1))

res_raw <- as.data.frame(results(dds, contrast = c("condition", "AD", "CN")))
res_shr <- as.data.frame(lfcShrink(dds, coef = "condition_AD_vs_CN", type = "normal", quiet = TRUE))

summ <- function(df, lab) {
  s <- subset(df, pvalue < 0.05 & abs(log2FoldChange) > 0.8)
  cat(sprintf("%-28s DEGs %4d  (up %4d / down %4d)  padj<0.05 %d\n",
      lab, nrow(s), sum(s$log2FoldChange > 0), sum(s$log2FoldChange < 0),
      sum(df$padj < 0.05, na.rm = TRUE)))
  s
}
cat("\n=== DEG counts: raw vs shrunken LFC (threshold p<0.05,|lfc|>0.8) ===\n")
cat(sprintf("Paper target:                DEGs  195  (up  100 / down  95)\n"))
s_raw <- summ(res_raw, "gene-level raw LFC:")
s_shr <- summ(res_shr, "gene-level shrunken LFC:")

# Recovery of paper's named up/down genes (from summary.md)
paper_down <- c("HTRA1","MELTF","KCNQ3","IGF1","TSPYL5","INPP5F","RPL23","RPL27A",
                "RPS19","ECHDC3","CKB","TRDN","IGFBP3","ITGB7","RPLP1","RPS12",
                "H3C","H4C5","H4C13","H4-16")
paper_up   <- c("PHLDA1","BCL2A1","POU5F1","PRKCB","NR4A2","NR4A3","FOSL1","SLIT3",
                "USP2","PLPP3","LILRA6")
chk <- function(genes, want, df) {
  d <- df[intersect(genes, rownames(df)), ]
  ok <- sum(sign(d$log2FoldChange) == want)
  cat(sprintf("  %-5s genes: %d/%d present, %d/%d correct direction\n",
      ifelse(want>0,"UP","DOWN"), nrow(d), length(genes), ok, nrow(d)))
}
cat("\n=== Paper marker-gene recovery (shrunken, all tested genes) ===\n")
chk(paper_down, -1, res_shr); chk(paper_up, 1, res_shr)

write.csv(res_shr, "../results/repro_rnaseq_gene_shrunk_all.csv")
write.csv(s_shr[order(s_shr$pvalue), ], "../results/repro_rnaseq_gene_shrunk_DEGs.csv")
cat("\nDone.\n")
