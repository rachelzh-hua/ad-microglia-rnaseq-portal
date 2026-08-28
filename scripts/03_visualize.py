#!/usr/bin/env python3
"""GSE332551 Reanalysis — Visualization
Volcano plots, sncRNA class breakdown, heatmaps
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
import os

os.makedirs('../figures', exist_ok=True)

# ============================================================
# 1. RNA-seq Volcano Plot
# ============================================================
rna = pd.read_csv('../results/rnaseq_deseq2_all.csv', index_col=0)
rna = rna.dropna(subset=['pvalue', 'log2FoldChange'])

fig, ax = plt.subplots(figsize=(8, 6))
neg_log10p = -np.log10(rna['pvalue'].clip(lower=1e-50))

# Color: significant DEGs
sig_up = (rna['pvalue'] < 0.05) & (rna['log2FoldChange'] > 0.8)
sig_down = (rna['pvalue'] < 0.05) & (rna['log2FoldChange'] < -0.8)
ns = ~(sig_up | sig_down)

ax.scatter(rna.loc[ns, 'log2FoldChange'], neg_log10p[ns], c='lightgrey', s=4, alpha=0.5, label='NS')
ax.scatter(rna.loc[sig_up, 'log2FoldChange'], neg_log10p[sig_up], c='#e74c3c', s=8, alpha=0.7, label=f'Up in AD ({sig_up.sum()})')
ax.scatter(rna.loc[sig_down, 'log2FoldChange'], neg_log10p[sig_down], c='#3498db', s=8, alpha=0.7, label=f'Down in AD ({sig_down.sum()})')

ax.axhline(-np.log10(0.05), ls='--', c='grey', lw=0.8)
ax.axvline(0.8, ls='--', c='grey', lw=0.8)
ax.axvline(-0.8, ls='--', c='grey', lw=0.8)

# Label top genes from the paper
paper_genes = ['HLA-DRA', 'AIF1', 'HTRA1', 'MELTF', 'POU5F1', 'PHLDA1', 'BCL2A1', 'PRKCB',
               'VCP', 'NANS', 'ACTN1', 'NPC2', 'IGF1', 'KCNQ3', 'INPP5F', 'TSPYL5']
for g in paper_genes:
    if g in rna.index and not pd.isna(rna.loc[g, 'pvalue']):
        row = rna.loc[g]
        if row['pvalue'] < 0.05 and abs(row['log2FoldChange']) > 0.8:
            ax.annotate(g, (row['log2FoldChange'], -np.log10(max(row['pvalue'], 1e-50))),
                       fontsize=6, ha='center', va='bottom')

ax.set_xlabel('log2 Fold Change (AD/CN)')
ax.set_ylabel('-log10(p-value)')
ax.set_title('RNA-seq: Sporadic AD vs CN iPSC-derived Microglia\n(p<0.05, |log2FC|>0.8)')
ax.legend(fontsize=8, loc='upper right')
plt.tight_layout()
plt.savefig('../figures/fig1_rnaseq_volcano.png', dpi=200)
plt.close()
print("Saved: fig1_rnaseq_volcano.png")

# ============================================================
# 2. sncRNA Volcano Plot
# ============================================================
snc = pd.read_csv('../results/sncrnaseq_deseq2_all.csv', index_col=0)
snc = snc.dropna(subset=['pvalue', 'log2FoldChange'])

def classify(name):
    if 'tRF' in str(name): return 'tRF'
    if 'miR' in str(name) or 'let-' in str(name): return 'miRNA'
    if 'piR' in str(name): return 'piRNA'
    if 'SNORD' in str(name) or 'wgRna' in str(name): return 'snoRNA'
    return 'other'

snc['type'] = [classify(n) for n in snc.index]

fig, ax = plt.subplots(figsize=(8, 6))
neg_log10p = -np.log10(snc['pvalue'].clip(lower=1e-50))
sig_up = (snc['pvalue'] < 0.05) & (snc['log2FoldChange'] > 1)
sig_down = (snc['pvalue'] < 0.05) & (snc['log2FoldChange'] < -1)
ns = ~(sig_up | sig_down)

colors = {'tRF': '#e67e22', 'miRNA': '#2ecc71', 'piRNA': '#9b59b6', 'snoRNA': '#1abc9c', 'other': 'grey'}

ax.scatter(snc.loc[ns, 'log2FoldChange'], neg_log10p[ns], c='lightgrey', s=6, alpha=0.3, label='NS')
for t in ['tRF', 'miRNA', 'piRNA', 'snoRNA']:
    mask_up = sig_up & (snc['type'] == t)
    mask_dn = sig_down & (snc['type'] == t)
    mask = mask_up | mask_dn
    if mask.sum() > 0:
        ax.scatter(snc.loc[mask, 'log2FoldChange'], neg_log10p[mask],
                  c=colors[t], s=15, alpha=0.8, label=f'{t} ({mask.sum()})')

ax.axhline(-np.log10(0.05), ls='--', c='grey', lw=0.8)
ax.axvline(1, ls='--', c='grey', lw=0.8)
ax.axvline(-1, ls='--', c='grey', lw=0.8)
ax.set_xlabel('log2 Fold Change (AD/CN)')
ax.set_ylabel('-log10(p-value)')
ax.set_title('sncRNA-seq: Sporadic AD vs CN iPSC-derived Microglia\n(p<0.05, |log2FC|>1)')
ax.legend(fontsize=8, loc='upper left')
plt.tight_layout()
plt.savefig('../figures/fig2_sncrnaseq_volcano.png', dpi=200)
plt.close()
print("Saved: fig2_sncrnaseq_volcano.png")

# ============================================================
# 3. sncRNA class breakdown (pie)
# ============================================================
sig_snc = pd.read_csv('../results/sncrnaseq_deseq2_DEsncRNAs.csv', index_col=0)

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

# 3a. All DE sncRNAs by type
if 'type' in sig_snc.columns:
    type_counts = sig_snc['type'].value_counts()
else:
    sig_snc['type'] = [classify(n) for n in sig_snc.index]
    type_counts = sig_snc['type'].value_counts()

clr = [colors.get(t, 'grey') for t in type_counts.index]
axes[0].pie(type_counts, labels=[f'{t}\n(n={c})' for t, c in type_counts.items()],
           colors=clr, autopct='%1.0f%%', startangle=90)
axes[0].set_title('DE sncRNA Classes')

# 3b. Direction by type
up_dn = sig_snc.groupby('type')['log2FoldChange'].apply(lambda x: pd.Series({
    'Up in AD': (x > 0).sum(), 'Down in AD': (x < 0).sum()
})).unstack()
up_dn.plot.bar(ax=axes[1], color=['#e74c3c', '#3498db'])
axes[1].set_title('DE sncRNAs: Direction by Class')
axes[1].set_ylabel('Count')
axes[1].set_xlabel('')
axes[1].legend(fontsize=8)
axes[1].tick_params(axis='x', rotation=0)

plt.tight_layout()
plt.savefig('../figures/fig3_sncRNA_breakdown.png', dpi=200)
plt.close()
print("Saved: fig3_sncRNA_breakdown.png")

# ============================================================
# 4. Top DEG heatmap (RNA-seq, top 40 by fold change)
# ============================================================
degs = pd.read_csv('../results/rnaseq_deseq2_DEGs.csv', index_col=0)
norm = pd.read_csv('../results/rnaseq_normalized_counts.csv', index_col=0)

# Top 40 by absolute FC
top40 = degs.reindex(degs['log2FoldChange'].abs().sort_values(ascending=False).index[:40])
top40_counts = norm.loc[norm.index.isin(top40.index)]
# Z-score per gene
top40_z = top40_counts.subtract(top40_counts.mean(axis=1), axis=0).div(top40_counts.std(axis=1), axis=0)
top40_z = top40_z.dropna()

fig, ax = plt.subplots(figsize=(5, 10))
# Cluster rows
if len(top40_z) > 2:
    row_link = linkage(pdist(top40_z.values, 'correlation'), method='average')
    row_order = dendrogram(row_link, no_plot=True)['leaves']
    top40_z = top40_z.iloc[row_order]

im = ax.imshow(top40_z.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
ax.set_yticks(range(len(top40_z)))
ax.set_yticklabels(top40_z.index, fontsize=6)
ax.set_xticks(range(len(top40_z.columns)))
ax.set_xticklabels(top40_z.columns, fontsize=8, rotation=45, ha='right')
ax.set_title('Top 40 DEGs (by |log2FC|)\nZ-scored Normalized Counts')
plt.colorbar(im, ax=ax, shrink=0.5, label='Z-score')
plt.tight_layout()
plt.savefig('../figures/fig4_rnaseq_heatmap_top40.png', dpi=200)
plt.close()
print("Saved: fig4_rnaseq_heatmap_top40.png")

# ============================================================
# 5. GO-relevant gene categories summary
# ============================================================
# Annotate paper-highlighted categories
go_cats = {
    'Extracellular/EV': ['CD302', 'FCGRT', 'OSCAR', 'SIRPB1', 'LRP8', 'MFSD2A'],
    'Chromatin/Histone': ['H3C', 'H4C5', 'H4C13', 'H4-16', 'EPOP'],
    'Cytoskeleton': ['ACTN1', 'TJP1', 'NES'],
    'Ribosomal': ['RPL23', 'RPL27A', 'RPLP1', 'RPS12', 'RPS19'],
    'Immune/MHC': ['HLA-DRA', 'AIF1', 'BCL2A1', 'CTLA4', 'IL12B'],
    'Signaling': ['PRKCB', 'NR4A2', 'NR4A3', 'FOSL1', 'IGF1']
}

fig, ax = plt.subplots(figsize=(8, 5))
cat_data = []
for cat, genes in go_cats.items():
    for g in genes:
        if g in degs.index:
            cat_data.append({'Category': cat, 'Gene': g,
                           'log2FC': degs.loc[g, 'log2FoldChange'],
                           'pvalue': degs.loc[g, 'pvalue']})
cat_df = pd.DataFrame(cat_data)
if len(cat_df) > 0:
    cat_colors = {'Extracellular/EV': '#3498db', 'Chromatin/Histone': '#e74c3c',
                  'Cytoskeleton': '#2ecc71', 'Ribosomal': '#f39c12',
                  'Immune/MHC': '#9b59b6', 'Signaling': '#1abc9c'}
    for cat in cat_df['Category'].unique():
        sub = cat_df[cat_df['Category'] == cat]
        ax.barh(sub['Gene'], sub['log2FC'], color=cat_colors.get(cat, 'grey'), label=cat, height=0.7)
    ax.set_xlabel('log2 Fold Change (AD/CN)')
    ax.set_title('Key DEGs by Functional Category')
    ax.axvline(0, color='black', lw=0.5)
    handles = [Patch(facecolor=cat_colors[c], label=c) for c in cat_df['Category'].unique()]
    ax.legend(handles=handles, fontsize=7, loc='lower right')
    plt.tight_layout()
    plt.savefig('../figures/fig5_go_categories.png', dpi=200)
    plt.close()
    print("Saved: fig5_go_categories.png")

# ============================================================
# 6. Paper comparison: overlap check
# ============================================================
# Paper reported 195 DEGs (95 down, 100 up)
paper_n = 195
our_n = len(degs)
# Check overlap with paper's table
paper_genes_full = [
    'ABCB1','ABCG1','ACD','ACVRL1','ADGRE3','AKR1B10','AMIGO2','ANGPTL2','ANKAR',
    'ARHGAP32','ARHGAP8','ARMCX2','ARV1','BBOF1','BCL2A1','C8B','CA10','CAMK2N2',
    'CATSPER1','CAVIN3','CCDC28B','CD302','CHRNA6','CKB','CLEC4F','CNGB1','CNKSR2',
    'COMT','CPB1','CPNE2','CREB3L4','CROCC','CTLA4','CYP3A7','DCST1','DERL3',
    'DIO1','DLX2','DRAM1','E2F5','ECHDC3','EDA','ENPP7','EPB41L1','EPOP','ERICH6',
    'ESRP2','FAM20A','FBXO47','FCGRT','FEZF2','FGF11','FLVCR2','FLYWCH2','FN3K',
    'FOSL1','FOXD4L1','GAD1','GIMAP5','GLCCI1','GLIS1','GNRHR','GPR35','GRID1',
    'GSTT1','GSTT2','GUCY2D','H3C','H4-16','H4C13','H4C5','HOXA1','HPDL','HSPA1A',
    'HTRA1','IER3','IFT140','IGF1','IGFBP3','IL12B','IL17C','INPP5F','ISL1','ITGB7',
    'JAG1','KBTBD13','KCNH8','KCNQ3','KRTAP9-2','LILRA6','LIN28B','LMOD2',
    'LPIN1','LRATD2','LRP8','LRRC4','LRRC43','MAEL','MAFF','MAN1A1','MARCHF9',
    'MCEE','MEI4','MELTF','MFSD2A','MND1','MUC2','MYLIP','NAGS','NEO1','NES',
    'NINJ2','NR4A2','NR4A3','OR10H4','OR2T8','OR2W3','OR4C13','OR8J1','OSCAR',
    'OSGIN1','PAM','PCDHGB2','PCDHGC3','PDCD5','PGAM2','PHETA2','PHLDA1','PHLDA2',
    'PLPP3','PMFBP1','POU5F1','POU6F2','PRAMEF12','PRKCB','PRR7','PYHIN1','RAB3IL1',
    'RANBP17','RASAL3','RBM41','RGS5','RIMBP3C','RIMS4','RNF145','RPL23','RPL27A',
    'RPLP1','RPS12','RPS19','RSPH10B2','RTN4RL1','RUNX3','SEMA4B','SEMG2','SERHL2',
    'SETSIP','SFN','SH3BGRL2','SIRPB1','SIRPD','SKOR1','SLC25A43','SLC39A11',
    'SLC43A1','SLIT3','SMARCD3','SMPDL3B','SNTG1','SNURF','SPATA31D','SPATA6L',
    'STARD8','SYCP2L','TBX6','TJP1','TMEM258','TMEM45A','TMSB15B','TRDN','TRIP6',
    'TSPYL5','TSTD1','TTC9','TUBA4A','TYMP','UBQLNL','UGT1A7','USP2','USP50',
    'ZFP3','ZNF728'
]

our_degs_set = set(degs.index)
paper_set = set(paper_genes_full)
overlap = our_degs_set & paper_set
only_paper = paper_set - our_degs_set
only_ours = our_degs_set - paper_set

summary_txt = f"""=== Reanalysis vs Paper Comparison ===
Paper DEGs: {len(paper_set)}
Our DEGs (p<0.05, |log2FC|>0.8): {len(our_degs_set)}
Overlap: {len(overlap)} ({100*len(overlap)/len(paper_set):.1f}% of paper DEGs recovered)
Paper-only: {len(only_paper)}
Reanalysis-only: {len(only_ours)}

sncRNA-seq:
Paper DE sncRNAs: 64
Our DE sncRNAs (p<0.05, |log2FC|>1): 68
"""
print(summary_txt)

with open('../results/comparison_summary.txt', 'w') as f:
    f.write(summary_txt)

print("All visualizations complete.")
