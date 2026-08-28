#!/usr/bin/env python3
"""GSE332551 — microglial state scoring (ssGSEA) on the corrected samples.

Rather than test single genes at n=5, score each sample against defined
microglial programs: homeostatic vs disease-associated (DAM/MGnD) vs
lipid-droplet-accumulating (LDAM), plus the SREBP cholesterol program, the
NF-kB inflammatory program, translation, and interferon. This turns the
GSEA themes into an interpretable, per-sample state read-out.

Writes results/microglia_state_scores.csv and figures/fig5_microglia_states.png.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
FIG = ROOT / "figures"

SIGS = {
    "Homeostatic": ["P2RY12", "TMEM119", "CX3CR1", "SALL1", "SELPLG", "CSF1R", "GPR34",
                    "MEF2C", "MAF", "SLC2A5", "SIGLEC11", "OLFML3"],
    "DAM / MGnD": ["TREM2", "APOE", "ITGAX", "CST7", "LPL", "CLEC7A", "SPP1", "GPNMB", "CD9",
                   "TYROBP", "AXL", "LGALS3", "CTSB", "CTSD", "FTH1", "CCL3", "CCL4", "CD63", "LILRB4"],
    "LDAM (lipid-droplet)": ["ACSL1", "DPYD", "NAMPT", "NCEH1", "NPC2", "SOAT1", "PLIN2",
                             "ABCA1", "GRN", "CD9"],
    "Cholesterol / SREBP": ["HMGCR", "HMGCS1", "SREBF2", "LDLR", "DHCR7", "DHCR24", "SQLE", "MVD",
                            "FDPS", "INSIG1", "MSMO1", "IDI1", "FDFT1", "CYP51A1", "LSS", "ACAT2", "MVK"],
    "Inflammatory / NF-kB": ["TNF", "IL1B", "NFKB1", "NFKBIA", "CXCL8", "CCL2", "IL6", "TNFAIP3",
                             "RELB", "BIRC3", "CXCL2", "PTGS2", "IER3", "SOD2"],
    "Translation / ribosome": ["RPL13", "RPL23", "RPS19", "RPS12", "RPLP1", "RPL27A", "EEF1A1",
                               "EEF2", "RPS6", "RPL7", "RPS3", "RPL11"],
    "Interferon": ["IFIT1", "IFIT3", "ISG15", "IRF7", "OAS1", "MX1", "STAT1", "IFITM3", "RSAD2", "USP18"],
}


def main() -> int:
    FIG.mkdir(exist_ok=True)
    expr = pd.read_csv(RES / "corrected_rnaseq_normalized_counts.csv", index_col=0)
    expr = np.log2(expr + 1)
    samples = list(expr.columns)
    groups = {s: ("AD" if s.endswith("_AD") else "CN") for s in samples}
    print(f"Samples: {samples}")

    ss = gp.ssgsea(data=expr, gene_sets=SIGS, sample_norm_method="rank",
                   min_size=5, outdir=None, threads=4)
    nes = ss.res2d.pivot(index="Term", columns="Name", values="NES").astype(float)
    nes = nes.reindex(index=list(SIGS), columns=samples)
    nes.to_csv(RES / "microglia_state_scores.csv")

    ad = [s for s in samples if groups[s] == "AD"]
    cn = [s for s in samples if groups[s] == "CN"]
    diff = (nes[ad].mean(axis=1) - nes[cn].mean(axis=1)).sort_values(ascending=False)
    print("\nAD - CN mean ssGSEA score (positive = higher in AD):")
    for k, v in diff.items():
        print(f"  {v:+.3f}  {k}")

    make_figure(nes, samples, groups, diff)
    return 0


def make_figure(nes, samples, groups, diff):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    # z-score each signature across samples for the heatmap
    z = nes.sub(nes.mean(axis=1), axis=0).div(nes.std(axis=1) + 1e-9, axis=0)
    order = diff.index.tolist()
    z = z.reindex(order)
    col_order = [s for s in samples if groups[s] == "AD"] + [s for s in samples if groups[s] == "CN"]
    z = z[col_order]

    fig, (axh, axb) = plt.subplots(1, 2, figsize=(11.5, 5.2), gridspec_kw={"width_ratios": [1.5, 1]})
    im = axh.imshow(z.values, cmap="RdBu_r", norm=TwoSlopeNorm(0, -2, 2), aspect="auto")
    axh.set_xticks(range(len(col_order)))
    axh.set_xticklabels([f"{s}" for s in col_order], rotation=45, ha="right", fontsize=9)
    axh.set_yticks(range(len(order))); axh.set_yticklabels(order, fontsize=9)
    for i, s in enumerate(col_order):
        axh.get_xticklabels()[i].set_color("#d64525" if groups[s] == "AD" else "#2f6db5")
    axh.set_title("Per-sample state score (row z-scored ssGSEA)", fontsize=10.5)
    fig.colorbar(im, ax=axh, fraction=0.046, pad=0.04, label="z-score")

    colors = ["#d64525" if v > 0 else "#2f6db5" for v in diff]
    axb.barh(range(len(diff)), diff.values, color=colors, edgecolor="white")
    axb.set_yticks(range(len(diff))); axb.set_yticklabels(diff.index, fontsize=9)
    axb.invert_yaxis(); axb.axvline(0, color="#7a8296", lw=0.8)
    axb.set_xlabel("AD − CN mean score")
    axb.set_title("Shift in AD (red = up in AD)", fontsize=10.5)
    for s in ("top", "right"):
        axb.spines[s].set_visible(False)
    fig.suptitle("GSE332551 — microglial state scoring (n = 2 AD, 3 CN; R17_AD excluded)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIG / "fig5_microglia_states.png", dpi=200)
    print("Wrote figures/fig5_microglia_states.png")


if __name__ == "__main__":
    sys.exit(main())
