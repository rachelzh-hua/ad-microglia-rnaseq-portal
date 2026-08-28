#!/usr/bin/env Rscript
# GSE332551 Reanalysis — RNA-seq DESeq2
# Sporadic AD vs CN iPSC-derived microglia

library(DESeq2)

# --- Load data ---
counts_raw <- read.csv("../data/GSE332551_counts_RNAseq.csv", row.names = 1)

# Collapse transcript-level to gene-level (sum isoforms)
gene_names <- sub("^[^:]+:", "", rownames(counts_raw))
counts_gene <- aggregate(counts_raw, by = list(Gene = gene_names), FUN = sum)
rownames(counts_gene) <- counts_gene$Gene
counts_gene$Gene <- NULL

# Sample metadata
coldata <- data.frame(
  sample = colnames(counts_gene),
  condition = ifelse(grepl("AD", colnames(counts_gene)), "AD", "CN"),
  row.names = colnames(counts_gene)
)
coldata$condition <- factor(coldata$condition, levels = c("CN", "AD"))

cat("Samples:\n")
print(coldata)
cat(sprintf("\nGenes (after collapsing isoforms): %d\n", nrow(counts_gene)))

# --- DESeq2 ---
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_gene),
  colData = coldata,
  design = ~ condition
)

# Pre-filter: keep genes with mean count >= 10
keep <- rowMeans(counts(dds)) >= 10
dds <- dds[keep, ]
cat(sprintf("Genes after filtering (mean >= 10): %d\n", sum(keep)))

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "AD", "CN"))

# --- Summarize ---
cat("\n=== DESeq2 Results Summary ===\n")
summary(res)

# Apply paper's thresholds: p < 0.05, |log2FC| > 0.8
sig <- subset(as.data.frame(res), pvalue < 0.05 & abs(log2FoldChange) > 0.8)
sig <- sig[order(sig$pvalue), ]
cat(sprintf("\nDEGs (p<0.05, |log2FC|>0.8): %d\n", nrow(sig)))
cat(sprintf("  Upregulated in AD: %d\n", sum(sig$log2FoldChange > 0)))
cat(sprintf("  Downregulated in AD: %d\n", sum(sig$log2FoldChange < 0)))

# Save results
write.csv(as.data.frame(res), "../results/rnaseq_deseq2_all.csv")
write.csv(sig, "../results/rnaseq_deseq2_DEGs.csv")

# Save normalized counts
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, "../results/rnaseq_normalized_counts.csv")

# Save size factors
cat("\nSize factors:\n")
print(sizeFactors(dds))

cat("\nDone. Results saved to ../results/\n")
