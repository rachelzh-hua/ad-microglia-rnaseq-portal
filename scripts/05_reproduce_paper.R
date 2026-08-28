#!/usr/bin/env Rscript
# GSE332551 — Reproduce the paper's RNA-seq DEG analysis
# Runs TRANSCRIPT-level (paper's approach) and GENE-level (original reanalysis)
# side by side so the 827-vs-195 discrepancy can be diagnosed.
# Paper: p < 0.05, |log2FC| > 0.8; reports 195 DEGs (100 up, 95 down).

suppressMessages(library(DESeq2))

counts_raw <- read.csv("../data/GSE332551_counts_RNAseq.csv", row.names = 1)
cat(sprintf("Raw transcripts: %d\n", nrow(counts_raw)))

coldata <- data.frame(
  condition = factor(ifelse(grepl("AD", colnames(counts_raw)), "AD", "CN"),
                     levels = c("CN", "AD")),
  row.names = colnames(counts_raw)
)

run_deseq <- function(mat, min_mean, label) {
  dds <- DESeqDataSetFromMatrix(round(mat), coldata, design = ~ condition)
  keep <- rowMeans(counts(dds)) >= min_mean
  dds <- dds[keep, ]
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("condition", "AD", "CN")))
  sig <- subset(res, pvalue < 0.05 & abs(log2FoldChange) > 0.8)
  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("Features after filter (mean>=%d): %d\n", min_mean, sum(keep)))
  cat(sprintf("DEG features (p<0.05,|lfc|>0.8): %d  (up %d / down %d)\n",
              nrow(sig), sum(sig$log2FoldChange > 0), sum(sig$log2FoldChange < 0)))
  cat(sprintf("  padj<0.05: %d\n", sum(res$padj < 0.05, na.rm = TRUE)))
  list(res = res, sig = sig, dds = dds)
}

# ---- Approach A: TRANSCRIPT level (paper) ----
tx <- run_deseq(counts_raw, 10, "TRANSCRIPT-level (paper approach)")
# Map significant transcripts -> unique gene symbols
tx_genes <- unique(sub("^[^:]+:", "", rownames(tx$sig)))
tx_sig <- tx$sig
tx_sig$gene <- sub("^[^:]+:", "", rownames(tx_sig))
# per-gene direction: take the most significant transcript per gene
tx_sig <- tx_sig[order(tx_sig$pvalue), ]
tx_gene_lvl <- tx_sig[!duplicated(tx_sig$gene), ]
cat(sprintf("  -> unique genes among DEG transcripts: %d  (up %d / down %d)\n",
            length(tx_genes),
            sum(tx_gene_lvl$log2FoldChange > 0),
            sum(tx_gene_lvl$log2FoldChange < 0)))

# ---- Approach B: GENE level (original reanalysis) ----
gene_names <- sub("^[^:]+:", "", rownames(counts_raw))
counts_gene <- aggregate(counts_raw, by = list(Gene = gene_names), FUN = sum)
rownames(counts_gene) <- counts_gene$Gene; counts_gene$Gene <- NULL
gn <- run_deseq(counts_gene, 10, "GENE-level (isoform-summed, original)")

# ---- Save transcript-level (paper-style) outputs ----
write.csv(tx$res, "../results/repro_rnaseq_transcript_all.csv")
write.csv(tx$sig, "../results/repro_rnaseq_transcript_DEGs.csv")
write.csv(tx_gene_lvl, "../results/repro_rnaseq_transcript_DEGs_byGene.csv")

cat("\nDone.\n")
