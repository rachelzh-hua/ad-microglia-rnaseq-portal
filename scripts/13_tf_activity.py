#!/usr/bin/env python3
"""GSE332551 — transcription-factor activity inference (decoupler + CollecTRI).

GSEA says *what* changed; this asks *who drives it*. We infer TF activity for
the AD-vs-CN contrast by running a univariate linear model (decoupler ULM) of
the gene-level Wald statistic against the signed CollecTRI regulons. Positive
activity = the TF's targets are coordinately up in AD (TF activated).

Writes results/tf_activity.csv and figures/fig6_tf_activity.png.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
import decoupler as dc

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
FIG = ROOT / "figures"

HIGHLIGHT = {"NFKB1", "RELA", "RELB", "FOS", "JUN", "JUNB", "FOSL1", "FOSL2",
             "SREBF1", "SREBF2", "STAT1", "NR1H3", "NR1H2", "MYC", "SPI1", "IRF7", "IRF8", "CEBPB"}


def main() -> int:
    FIG.mkdir(exist_ok=True)
    df = pd.read_csv(RES / "corrected_rnaseq_deseq2_all.csv").dropna(subset=["log2FC_raw", "lfcSE"])
    df = df[df["lfcSE"] > 0]
    df["stat"] = df["log2FC_raw"] / df["lfcSE"]
    df = df.sort_values("stat", key=lambda s: s.abs(), ascending=False).drop_duplicates("gene")
    mat = pd.DataFrame([df.set_index("gene")["stat"]], index=["AD_vs_CN"])

    net = dc.op.collectri(organism="human")
    acts, pvals = dc.mt.ulm(data=mat, net=net, tmin=5)

    out = pd.DataFrame({
        "TF": acts.columns,
        "activity": acts.loc["AD_vs_CN"].values,
        "pvalue": pvals.loc["AD_vs_CN"].values,
    })
    out["FDR"] = false_discovery_control(out["pvalue"].fillna(1.0), method="bh")
    out["direction"] = np.where(out["activity"] >= 0, "activated_in_AD", "repressed_in_AD")
    out = out.sort_values("activity", ascending=False).reset_index(drop=True)
    out.to_csv(RES / "tf_activity.csv", index=False)

    sig = out[out["FDR"] < 0.05]
    print(f"{out.shape[0]} TFs scored; {len(sig)} at FDR<0.05")
    print("\nTop ACTIVATED in AD:")
    for _, r in out.head(10).iterrows():
        star = " *" if r.FDR < 0.05 else ""
        print(f"  {r.activity:+.2f}  FDR {r.FDR:.1e}  {r.TF}{star}")
    print("Top REPRESSED in AD:")
    for _, r in out.tail(8).iloc[::-1].iterrows():
        star = " *" if r.FDR < 0.05 else ""
        print(f"  {r.activity:+.2f}  FDR {r.FDR:.1e}  {r.TF}{star}")

    make_figure(out)
    return 0


def make_figure(out: pd.DataFrame):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    UP, DOWN = "#d64525", "#2f6db5"
    top = pd.concat([out.head(13), out.tail(13)]).drop_duplicates("TF")
    top = top.sort_values("activity")
    colors = [UP if v > 0 else DOWN for v in top["activity"]]

    fig, ax = plt.subplots(figsize=(8.4, 8.2))
    ax.barh(range(len(top)), top["activity"], color=colors, edgecolor="white")
    ax.set_yticks(range(len(top)))
    labels = []
    for _, r in top.iterrows():
        mark = "  ★" if r.FDR < 0.05 else ""
        labels.append(f"{r.TF}{mark}")
    ax.set_yticklabels(labels, fontsize=9)
    for i, (_, r) in enumerate(top.iterrows()):
        if r.TF in HIGHLIGHT:
            ax.get_yticklabels()[i].set_fontweight("bold")
    ax.axvline(0, color="#7a8296", lw=0.8)
    ax.set_xlabel("Inferred TF activity (AD vs CN)")
    ax.set_title("GSE332551 — TF-activity inference (decoupler ULM · CollecTRI)\n"
                 "red = activated in AD · blue = repressed · ★ FDR < 0.05 · bold = expected regulator",
                 fontsize=10.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "fig6_tf_activity.png", dpi=200)
    print("Wrote figures/fig6_tf_activity.png")


if __name__ == "__main__":
    sys.exit(main())
