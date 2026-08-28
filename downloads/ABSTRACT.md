# GSE332551 Reanalysis: Project Abstract

**An outlier-corrected reanalysis of RNA-seq and small-RNA-seq in iPSC-derived microglia from sporadic Alzheimer's disease**

Re-analysis of Wu, Choi, Liu et al., *Frontiers in Neuroscience* (2026), GEO **GSE332551**

---

**Background.** Alzheimer's disease (AD) reshapes microglial gene and small-RNA programs, but patient-derived models are scarce. The original study profiled iPSC-derived microglia-like cells (iMG) from 3 sporadic AD donors and 3 cognitively normal (CN) controls, matched for age and sex and all APOE3/3, using total RNA-seq and T4 PNK small-RNA-seq. We reanalyzed the deposited count matrices to check whether the reported differential expression reproduces from public data and to put the DEG calls on a firmer footing.

**Methods.** We ran DESeq2 (v1.46) on both matrices (AD vs CN, n = 3 per group). RNA-seq transcript counts were collapsed to genes, leaving 13,611 genes after a mean-count filter of 10. Small-RNA features were tested directly (1,467 kept at mean ≥ 5) and classified as piRNA, tRF, miRNA, or snoRNA. We added normalization checks (library size, DESeq2 size factors, inter-sample correlation, PCA), applied log-fold-change shrinkage with FDR control, and re-ran the RNA-seq analysis with the outlier removed.

**Results.** The small-RNA results reproduced well. We found 68 differential sncRNAs, close to the 64 the paper reported, with a fold-change correlation of r = 1.000 and the same class pattern: piRNAs (35) and tRFs (21) up in AD, and miRNAs and snoRNAs mostly down. The RNA-seq results did not reproduce as cleanly. The deposited matrix gave 827 DEGs with a 3.6-to-1 skew toward "up in AD", and that skew traced back to one sample, R17_AD. It had 3.7 M reads against 8.4 to 15.2 M for the others, a size factor of 0.44, a genome-wide Spearman correlation near 0.80 where the others sat around 0.95, and a lone position on PCA. Dropping R17_AD removed the skew, leaving 387 DEGs split 134 up and 253 down. It also raised the FDR-significant genes from 11 to 30 and left 13 high-confidence DEGs by FDR and shrunken fold-change, among them PHLDA1 up, SLC44A2 down, ALK up, and KCNQ3 down. Every named marker gene from the paper came back in the right direction (11 of 11 up, 18 of 18 down).

**Conclusions.** The underlying biology holds up. The small-RNA class shifts and the marker-gene directions both reproduce. The exact RNA-seq DEG count does not, because it is unstable at n = 3 and was inflated by one low-depth sample. We recommend working from an outlier-corrected, FDR-based DEG list. The results, figures, analysis scripts, and an interactive DEG table that you can filter by threshold are all provided in this portal.

*Prepared Reanalysis completed 2026-07.*
