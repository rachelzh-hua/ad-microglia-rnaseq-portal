#!/usr/bin/env Rscript
# GSE332551 — Corrected RNA-seq DEG list
# Definitive reanalysis: gene-level DESeq2, R17_AD outlier excluded,
# normal LFC shrinkage, FDR reported. Supersedes the original 827-DEG list.
suppressMessages(library(DESeq2))

counts_raw <- read.csv("../data/GSE332551_counts_RNAseq.csv", row.names = 1)
gene_names <- sub("^[^:]+:", "", rownames(counts_raw))
counts_gene <- aggregate(counts_raw, by = list(Gene = gene_names), FUN = sum)
rownames(counts_gene) <- counts_gene$Gene; counts_gene$Gene <- NULL

# --- Exclude the R17_AD outlier (see REPRODUCTION_REPORT.md) ---
drop <- "R17_AD"
counts_gene <- counts_gene[, !colnames(counts_gene) %in% drop]
cat(sprintf("Samples used: %s\n", paste(colnames(counts_gene), collapse = ", ")))

coldata <- data.frame(
  condition = factor(ifelse(grepl("AD", colnames(counts_gene)), "AD", "CN"),
                     levels = c("CN", "AD")),
  row.names = colnames(counts_gene))

dds <- DESeqDataSetFromMatrix(round(counts_gene), coldata, design = ~ condition)
dds <- dds[rowMeans(counts(dds)) >= 10, ]
cat(sprintf("Genes after filter (mean>=10): %d\n", nrow(dds)))
dds <- DESeq(dds, quiet = TRUE)

# Raw (unshrunken) results = paper-comparable magnitudes + p/padj
raw <- as.data.frame(results(dds, contrast = c("condition", "AD", "CN")))
# Shrunken fold changes = stable magnitudes, robust to low counts
shr <- as.data.frame(lfcShrink(dds, coef = "condition_AD_vs_CN",
                               type = "normal", quiet = TRUE))

res <- data.frame(
  gene           = rownames(raw),
  baseMean       = raw$baseMean,
  log2FC_raw     = raw$log2FoldChange,
  log2FC_shrunk  = shr$log2FoldChange[match(rownames(raw), rownames(shr))],
  lfcSE          = raw$lfcSE,
  pvalue         = raw$pvalue,
  padj           = raw$padj,
  row.names      = NULL, stringsAsFactors = FALSE)
res$direction <- ifelse(res$log2FC_raw > 0, "up_in_AD", "down_in_AD")

# --- Two DEG tiers ---
# Tier 1 (paper-comparable, exploratory): raw p<0.05 & |raw log2FC|>0.8
deg_expl <- subset(res, !is.na(pvalue) & pvalue < 0.05 & abs(log2FC_raw) > 0.8)
deg_expl <- deg_expl[order(deg_expl$pvalue), ]
# Tier 2 (high-confidence): FDR<0.05 & |shrunken log2FC|>0.8
deg_hc <- subset(res, !is.na(padj) & padj < 0.05 & abs(log2FC_shrunk) > 0.8)
deg_hc <- deg_hc[order(deg_hc$padj), ]

rpt <- function(s, lab, fc)
  cat(sprintf("%-40s %4d  (up %3d / down %3d)\n", lab, nrow(s),
      sum(s[[fc]] > 0), sum(s[[fc]] < 0)))
cat("\n=== Corrected DEG counts (R17_AD excluded) ===\n")
cat(sprintf("%-40s %4s  (up %3s / down %3s)\n", "Paper (for reference):", "195", "100", " 95"))
rpt(deg_expl, "Tier 1 exploratory (raw p<0.05,|FC|>0.8):", "log2FC_raw")
rpt(deg_hc,   "Tier 2 high-conf (FDR<0.05,|shrunkFC|>0.8):", "log2FC_shrunk")

write.csv(res,      "../results/corrected_rnaseq_deseq2_all.csv",    row.names = FALSE)
write.csv(deg_expl, "../results/corrected_rnaseq_DEGs_exploratory.csv", row.names = FALSE)
write.csv(deg_hc,   "../results/corrected_rnaseq_DEGs_highconf.csv",  row.names = FALSE)

cat("\nTier 2 high-confidence DEGs (FDR<0.05):\n")
print(deg_hc[, c("gene","baseMean","log2FC_raw","log2FC_shrunk","padj","direction")],
      row.names = FALSE, digits = 3)
cat("\nSaved corrected_rnaseq_* to results/\n")
