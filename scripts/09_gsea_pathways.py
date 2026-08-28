#!/usr/bin/env python3
"""GSE332551 — preranked GSEA on the corrected RNA-seq analysis.

Threshold-free pathway analysis. Instead of testing an (unstable, n=3) DEG
list, rank every gene by its Wald statistic (log2FC / lfcSE) and run GSEA
against MSigDB Hallmark, GO-BP, Reactome and KEGG. Enrichments that survive
here do not depend on any p-value cutoff, so they sidestep the R17_AD /
threshold instability entirely.

Writes results/gsea_pathways.csv and figures/fig3_gsea.png.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
FIG = ROOT / "figures"

LIBRARIES = {
    "Hallmark": "MSigDB_Hallmark_2020",
    "GO-BP": "GO_Biological_Process_2023",
    "Reactome": "Reactome_2022",
    "KEGG": "KEGG_2021_Human",
}


def build_ranking() -> pd.DataFrame:
    df = pd.read_csv(RES / "corrected_rnaseq_deseq2_all.csv")
    df = df.dropna(subset=["log2FC_raw", "lfcSE"])
    df = df[df["lfcSE"] > 0]
    df["stat"] = df["log2FC_raw"] / df["lfcSE"]          # Wald statistic
    df = df.dropna(subset=["stat"])
    # collapse duplicate symbols to the most extreme statistic
    df["absstat"] = df["stat"].abs()
    df = df.sort_values("absstat", ascending=False).drop_duplicates("gene")
    rnk = df[["gene", "stat"]].sort_values("stat", ascending=False).reset_index(drop=True)
    return rnk


def main() -> int:
    FIG.mkdir(exist_ok=True)
    rnk = build_ranking()
    print(f"Ranked {len(rnk)} genes (stat range {rnk['stat'].min():.2f} … {rnk['stat'].max():.2f})")

    frames = []
    for label, lib in LIBRARIES.items():
        print(f"GSEA vs {label} ({lib}) …")
        try:
            pre = gp.prerank(
                rnk=rnk, gene_sets=lib, min_size=10, max_size=500,
                permutation_num=1000, seed=42, threads=4, no_plot=True, verbose=False,
            )
        except Exception as e:  # noqa: BLE001
            print(f"  skipped {label}: {e}")
            continue
        res = pre.res2d.copy()
        res["Library"] = label
        frames.append(res)

    if not frames:
        print("No GSEA results (Enrichr unreachable?).")
        return 1

    out = pd.concat(frames, ignore_index=True)
    for col in ("NES", "NOM p-val", "FDR q-val"):
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out = out.rename(columns={"Term": "pathway"})
    out = out[["Library", "pathway", "ES", "NES", "NOM p-val", "FDR q-val", "Lead_genes"]]
    out = out.sort_values("FDR q-val").reset_index(drop=True)
    out.to_csv(RES / "gsea_pathways.csv", index=False)
    sig = out[out["FDR q-val"] < 0.25]
    print(f"Wrote gsea_pathways.csv: {len(out)} sets tested, {len(sig)} at FDR<0.25")

    make_figure(sig if len(sig) else out)
    return 0


def make_figure(df: pd.DataFrame):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    INDIGO, UP, DOWN = "#512feb", "#d64525", "#2f6db5"
    top = df.reindex(df["NES"].abs().sort_values(ascending=False).index).head(16)
    top = top.sort_values("NES")
    labels = [f"{p[:46]}  ·  {lib}" for p, lib in zip(top["pathway"], top["Library"])]
    colors = [UP if v > 0 else DOWN for v in top["NES"]]

    fig, ax = plt.subplots(figsize=(9.6, max(4, 0.42 * len(top) + 1.3)))
    ax.barh(range(len(top)), top["NES"], color=colors, edgecolor="white")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.axvline(0, color="#7a8296", lw=0.8)
    ax.set_xlabel("Normalized enrichment score (AD vs CN)")
    ax.set_title("GSE332551 corrected RNA-seq — GSEA (FDR < 0.25)\nred = up in AD, blue = down in AD",
                 fontsize=11)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "fig3_gsea.png", dpi=200)
    print("Wrote figures/fig3_gsea.png")


if __name__ == "__main__":
    sys.exit(main())
