#!/usr/bin/env python3
"""GSE332551 — Figures from the CORRECTED analysis (R17_AD excluded, shrunken LFC, FDR)."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

RES = "../results/"
FIG = "../figures/"

res = pd.read_csv(RES + "corrected_rnaseq_deseq2_all.csv")
hc  = pd.read_csv(RES + "corrected_rnaseq_DEGs_highconf.csv")
norm = pd.read_csv(RES + "corrected_rnaseq_normalized_counts.csv", index_col=0)

FC, P = "log2FC_shrunk", "padj"
UP, DOWN, NS, HL = "#c0392b", "#2471a3", "#b0b0b0", "#000000"

# ============================================================
# FIG 1 — Volcano (shrunken log2FC vs FDR), high-conf genes labeled
# ============================================================
d = res.dropna(subset=[FC, P]).copy()
d["nlp"] = -np.log10(d[P].clip(lower=1e-300))
sig = (d[P] < 0.05) & (d[FC].abs() > 0.8)
up, dn = sig & (d[FC] > 0), sig & (d[FC] < 0)

fig, ax = plt.subplots(figsize=(8, 6.5))
ax.scatter(d.loc[~sig, FC], d.loc[~sig, "nlp"], s=6, c=NS, alpha=0.35, lw=0)
ax.scatter(d.loc[up, FC], d.loc[up, "nlp"], s=10, c=UP, alpha=0.8, lw=0, label=f"Up in AD ({int(up.sum())})")
ax.scatter(d.loc[dn, FC], d.loc[dn, "nlp"], s=10, c=DOWN, alpha=0.8, lw=0, label=f"Down in AD ({int(dn.sum())})")
ax.axhline(-np.log10(0.05), color="k", ls="--", lw=0.6, alpha=0.6)
ax.axvline(0.8, color="k", ls=":", lw=0.5, alpha=0.5); ax.axvline(-0.8, color="k", ls=":", lw=0.5, alpha=0.5)

paper_markers = {"PHLDA1","KCNQ3","CKB","LRP8","SLIT3","LRRC4"}
hcd = d[d["gene"].isin(hc["gene"])].copy()
ymax = d["nlp"].max()
# Stack labels in a vertical column on each side with leader lines -> no overlap
def place(side):
    sub = hcd[(hcd[FC] > 0) if side == "right" else (hcd[FC] < 0)].sort_values("nlp", ascending=False)
    n = len(sub)
    tx = d[FC].max() * 1.02 if side == "right" else d[FC].min() * 1.02
    ha = "left" if side == "right" else "right"
    ys = np.linspace(ymax * 0.98, ymax * 0.30, n) if n > 1 else [ymax * 0.6]
    for (_, r), ty in zip(sub.iterrows(), ys):
        ax.annotate(r["gene"], xy=(r[FC], r["nlp"]), xytext=(tx, ty),
                    fontsize=7, va="center", ha=ha,
                    fontweight="bold" if r["gene"] in paper_markers else "normal",
                    color=HL if r["gene"] in paper_markers else "#444",
                    arrowprops=dict(arrowstyle="-", color="#999", lw=0.4,
                                    shrinkA=0, shrinkB=2))
place("right"); place("left")
ax.set_xlim(d[FC].min() * 1.25, d[FC].max() * 1.25)
ax.set_xlabel("log2 fold change (shrunken), AD vs CN")
ax.set_ylabel("-log10(FDR)")
ax.set_title("GSE332551 corrected RNA-seq volcano\n(R17_AD outlier excluded; bold = confirmed in paper)")
ax.legend(loc="upper center", fontsize=9, frameon=False)
plt.tight_layout(); plt.savefig(FIG + "fig1_corrected_volcano.png", dpi=200); plt.close()
print("saved fig1_corrected_volcano.png")

# ============================================================
# FIG 2 — Heatmap of the 13 high-confidence DEGs (z-scored norm counts)
# ============================================================
genes = hc.sort_values(FC, ascending=False)["gene"].tolist()
mat = norm.loc[[g for g in genes if g in norm.index]]
z = mat.apply(lambda r: (r - r.mean()) / (r.std(ddof=0) + 1e-9), axis=1)
# order columns: AD then CN
cols = [c for c in z.columns if "AD" in c] + [c for c in z.columns if "CN" in c]
z = z[cols]

cmap = LinearSegmentedColormap.from_list("bwr2", [DOWN, "#f7f7f7", UP])
fig, ax = plt.subplots(figsize=(6, 6))
im = ax.imshow(z.values, aspect="auto", cmap=cmap, vmin=-1.5, vmax=1.5)
ax.set_yticks(range(len(z))); ax.set_yticklabels(z.index, fontsize=8)
ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=45, ha="right", fontsize=8)
for i, c in enumerate(cols):
    ax.get_xticklabels()[i].set_color(UP if "AD" in c else DOWN)
ax.set_title("13 high-confidence DEGs (FDR<0.05)\nrow z-scored normalized counts", fontsize=10)
cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04); cb.set_label("z-score", fontsize=8)
plt.tight_layout(); plt.savefig(FIG + "fig2_corrected_heatmap.png", dpi=200); plt.close()
print("saved fig2_corrected_heatmap.png")
