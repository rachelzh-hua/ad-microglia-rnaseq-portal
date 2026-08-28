#!/usr/bin/env python3
"""GSE332551 — miRNA -> mRNA integration across the two deposited assays.

If the differential miRNAs are functional, their mRNA targets should move the
opposite way in the RNA-seq: targets of DOWN-in-AD miRNAs should be de-repressed
(up), targets of UP-in-AD miRNAs should be repressed (down). We test this by
comparing the RNA-seq Wald statistic of each miRNA's validated targets
(miRTarBase) against the genome background (one-sided Mann-Whitney U).

Writes results/mirna_mrna_integration.csv and figures/fig4_mirna_mrna.png.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import gseapy as gp

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
FIG = ROOT / "figures"
MIN_TARGETS = 8


def ranking() -> dict:
    df = pd.read_csv(RES / "corrected_rnaseq_deseq2_all.csv").dropna(subset=["log2FC_raw", "lfcSE"])
    df = df[df["lfcSE"] > 0]
    df["stat"] = df["log2FC_raw"] / df["lfcSE"]
    df = df.sort_values("stat", key=lambda s: s.abs(), ascending=False).drop_duplicates("gene")
    return dict(zip(df["gene"], df["stat"]))


def target_lookup(lib: dict):
    """map many key spellings -> target set (handles arm suffix drops)."""
    norm = {}
    for k, v in lib.items():
        s = set(v)
        key = k.lower()
        base = key.rsplit("-", 1)[0]
        norm[key] = norm.get(key, set()) | s
        norm[base] = norm.get(base, set()) | s
    def get(mir):
        m = mir.lower()
        return norm.get(m) or norm.get(m.rsplit("-", 1)[0]) or set()
    return get


def main() -> int:
    FIG.mkdir(exist_ok=True)
    stat = ranking()
    genes = set(stat)
    bg = np.array(list(stat.values()))

    de = pd.read_csv(RES / "sncrnaseq_DE_all.csv")
    mir = de[(de["class"] == "miRNA") & (de.pvalue < 0.05) & (de.log2FoldChange.abs() > 1)].copy()

    lib = gp.get_library("miRTarBase_2017")
    get_targets = target_lookup(lib)

    rows = []
    pooled = {"down_in_AD": set(), "up_in_AD": set()}
    for _, r in mir.iterrows():
        tset = get_targets(r.sncRNA) & genes
        if len(tset) < MIN_TARGETS:
            rows.append(dict(miRNA=r.sncRNA, mir_dir=r.direction, mir_log2FC=round(r.log2FoldChange, 2),
                             n_targets=len(tset), target_median_stat=np.nan, expected="", concordant="",
                             mwu_p=np.nan))
            continue
        tvals = np.array([stat[g] for g in tset])
        # DOWN miRNA -> targets up (greater); UP miRNA -> targets down (less)
        alt = "greater" if r.direction == "down_in_AD" else "less"
        p = mannwhitneyu(tvals, bg, alternative=alt).pvalue
        med = float(np.median(tvals))
        concordant = (med > 0) if r.direction == "down_in_AD" else (med < 0)
        pooled[r.direction] |= tset
        rows.append(dict(miRNA=r.sncRNA, mir_dir=r.direction, mir_log2FC=round(r.log2FoldChange, 2),
                         n_targets=len(tset), target_median_stat=round(med, 3),
                         expected="targets up" if alt == "greater" else "targets down",
                         concordant=bool(concordant), mwu_p=p))
    out = pd.DataFrame(rows).sort_values(["mir_dir", "mwu_p"])
    out.to_csv(RES / "mirna_mrna_integration.csv", index=False)

    # pooled directional tests
    def pooled_test(tset, alt):
        t = np.array([stat[g] for g in tset])
        return len(tset), float(np.median(t)), mannwhitneyu(t, bg, alternative=alt).pvalue
    nd, md, pd_ = pooled_test(pooled["down_in_AD"], "greater")
    nu, mu, pu = pooled_test(pooled["up_in_AD"], "less")

    tested = out.dropna(subset=["mwu_p"])
    print(f"Tested {len(tested)} DE miRNAs with >= {MIN_TARGETS} targets.")
    print(f"  concordant direction: {tested['concordant'].sum()}/{len(tested)}")
    print(f"  DOWN-miRNA targets (n={nd}) median stat {md:+.3f}, up-shift p={pd_:.2e}")
    print(f"  UP-miRNA targets   (n={nu}) median stat {mu:+.3f}, down-shift p={pu:.2e}")

    make_figure(stat, pooled, bg, (nd, md, pd_), (nu, mu, pu))
    return 0


def make_figure(stat, pooled, bg, down, up):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    UP, DOWN, GREY = "#d64525", "#2f6db5", "#9aa1b4"
    dvals = [stat[g] for g in pooled["down_in_AD"]]
    uvals = [stat[g] for g in pooled["up_in_AD"]]
    groups = [("targets of\ndown-in-AD miRNAs\n(expect ↑)", dvals, DOWN),
              ("all genes\n(background)", list(bg), GREY),
              ("targets of\nup-in-AD miRNAs\n(expect ↓)", uvals, UP)]

    fig, ax = plt.subplots(figsize=(8, 5))
    parts = ax.violinplot([g[1] for g in groups], showmedians=True, widths=0.8)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(groups[i][2]); pc.set_alpha(0.55); pc.set_edgecolor("white")
    for key in ("cmedians", "cbars", "cmins", "cmaxes"):
        parts[key].set_color("#4a5163")
    ax.axhline(0, color="#7a8296", lw=0.8, ls="--")
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels([g[0] for g in groups], fontsize=9)
    ax.set_ylabel("RNA-seq Wald statistic (AD vs CN)")
    ax.set_ylim(-6, 8)
    ax.set_title(f"miRNA → mRNA integration\ndown-miRNA targets shift up (p={down[2]:.1e}), "
                 f"up-miRNA targets shift down (p={up[2]:.1e})", fontsize=10.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "fig4_mirna_mrna.png", dpi=200)
    print("Wrote figures/fig4_mirna_mrna.png")


if __name__ == "__main__":
    sys.exit(main())
