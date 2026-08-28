#!/usr/bin/env Rscript
# GSE332551 Reanalysis — sncRNA-seq DESeq2
# Sporadic AD vs CN iPSC-derived microglia

library(DESeq2)

# --- Load data ---
counts_raw <- read.csv("../data/GSE332551_counts_sncRNAseq.csv", row.names = 1)

# Sample metadata
coldata <- data.frame(
  sample = colnames(counts_raw),
  condition = ifelse(grepl("AD", colnames(counts_raw)), "AD", "CN"),
  row.names = colnames(counts_raw)
)
coldata$condition <- factor(coldata$condition, levels = c("CN", "AD"))

cat("Samples:\n")
print(coldata)
cat(sprintf("\nsncRNAs total: %d\n", nrow(counts_raw)))

# --- DESeq2 ---
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_raw),
  colData = coldata,
  design = ~ condition
)

# Pre-filter: keep sncRNAs with mean count >= 10 in either group
keep <- rowMeans(counts(dds)) >= 5  # slightly relaxed for small RNAs
dds <- dds[keep, ]
cat(sprintf("sncRNAs after filtering (mean >= 5): %d\n", sum(keep)))

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "AD", "CN"))

# --- Summarize ---
cat("\n=== DESeq2 Results Summary ===\n")
summary(res)

# Apply paper's thresholds: p < 0.05, |log2FC| > 1
sig <- subset(as.data.frame(res), pvalue < 0.05 & abs(log2FoldChange) > 1)
sig <- sig[order(sig$pvalue), ]
cat(sprintf("\nDE sncRNAs (p<0.05, |log2FC|>1): %d\n", nrow(sig)))

# Classify by type
classify_sncrna <- function(name) {
  if (grepl("^tRF", name)) return("tRF")
  if (grepl("^hsa-miR|^hsa-let", name)) return("miRNA")
  if (grepl("^hsa-piR", name)) return("piRNA")
  if (grepl("SNORD|wgRna", name)) return("snoRNA")
  return("other")
}
if (nrow(sig) > 0) {
  sig$type <- sapply(rownames(sig), classify_sncrna)
  cat("\nDE sncRNAs by type:\n")
  print(table(sig$type, sig$log2FoldChange > 0))
}

# Save results
write.csv(as.data.frame(res), "../results/sncrnaseq_deseq2_all.csv")
write.csv(sig, "../results/sncrnaseq_deseq2_DEsncRNAs.csv")

# Save normalized counts
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, "../results/sncrnaseq_normalized_counts.csv")

cat("\nDone. Results saved to ../results/\n")
