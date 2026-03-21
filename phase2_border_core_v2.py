#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 2-B — Border vs Core enrichment analysis (GSE144735)
Goal: Test whether Pre-driver myeloid candidates are enriched at the tumor border
      (supporting the "physical barrier" hypothesis)
Data: processed_data/GSE144735.h5ad  (already has TAM_score, Predriver_score)
Output: phase2_results/spatial_v2/  +  figures/phase2/spatial_v2/
"""
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import os, sys, io, warnings
warnings.filterwarnings('ignore')

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42
plt.rcParams['font.family']  = 'Arial'

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

RES = 'phase2_results/spatial_v2'
FIG = 'figures/phase2/spatial_v2'
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

print('=' * 65)
print('Phase 2-B  Border vs Core (GSE144735)')
print('=' * 65)

# ─────────────────────────────────────────────────────────────────
# 1. Load GSE144735
# ─────────────────────────────────────────────────────────────────
print('\n[1/5] Loading GSE144735 ...')
adata = sc.read_h5ad('processed_data/GSE144735.h5ad')
print(f'  Shape: {adata.shape}')
print(f'  Class distribution:\n{adata.obs["Class"].value_counts().to_string()}')
print(f'  Cell types:\n{adata.obs["Cell_type"].value_counts().to_string()}')

# Re-score with enriched Pre-driver gene set (SPP1+ TAM signature)
PREDRIVER_GENES = ['SPP1', 'APOE', 'C1QC', 'LGALS3', 'MRC1', 'VEGFA',
                   'FTH1', 'HMOX1', 'SOD2', 'CD68', 'HLA-DRA', 'HLA-DRB1',
                   'FCGR1A', 'CTSB', 'RGS1', 'GPR183', 'CD274', 'IL10',
                   'TGFB1', 'SOCS3']
INFLAM_GENES = ['S100A8', 'S100A9', 'CXCL8', 'IL1B', 'CCL3', 'CCL4',
                'TNF', 'IL6', 'NLRP3', 'G0S2', 'SELL', 'CD14', 'FCGR3A']

avail_pd = [g for g in PREDRIVER_GENES if g in adata.var_names]
avail_if = [g for g in INFLAM_GENES    if g in adata.var_names]
print(f'\n  Pre-driver genes available: {len(avail_pd)}/{len(PREDRIVER_GENES)}')
print(f'  Inflammatory genes available: {len(avail_if)}/{len(INFLAM_GENES)}')

sc.tl.score_genes(adata, avail_pd, score_name='Predriver_v2', use_raw=False)
sc.tl.score_genes(adata, avail_if, score_name='Inflam_score', use_raw=False)

# ─────────────────────────────────────────────────────────────────
# 2. Myeloid-only analysis
# ─────────────────────────────────────────────────────────────────
print('\n[2/5] Myeloid analysis ...')
my = adata[adata.obs['Cell_type'] == 'Myeloids'].copy()
print(f'  Myeloid cells: {my.n_obs:,}')
print(pd.crosstab(my.obs['Cell_subtype'], my.obs['Class']).to_string())

# ─────────────────────────────────────────────────────────────────
# 3. Statistics: Pre-driver score by Class
# ─────────────────────────────────────────────────────────────────
print('\n[3/5] Statistical tests ...')
CLASS_ORDER = ['Normal', 'Border', 'Tumor']
CLASS_COLORS = {'Normal': '#4DBBD5', 'Border': '#00A087', 'Tumor': '#E64B35'}

def kruskal_pairwise(data_dict):
    """Kruskal-Wallis + Mann-Whitney pairwise for 3 groups."""
    groups = list(data_dict.keys())
    vals   = list(data_dict.values())
    H, p_kw = stats.kruskal(*vals)
    pairs = {}
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            _, p = stats.mannwhitneyu(vals[i], vals[j], alternative='two-sided')
            pairs[f'{groups[i]}_vs_{groups[j]}'] = p
    return H, p_kw, pairs

results = {}
for score_col in ['TAM_score', 'CD8Tex_score', 'Predriver_score', 'Predriver_v2', 'Inflam_score']:
    for cell_grp, ad in [('All', adata), ('Myeloid', my)]:
        data = {cls: ad.obs.loc[ad.obs['Class'] == cls, score_col].values
                for cls in CLASS_ORDER}
        H, p_kw, pairs = kruskal_pairwise(data)
        for cls in CLASS_ORDER:
            results[(score_col, cell_grp, cls)] = {
                'mean': data[cls].mean(),
                'median': np.median(data[cls]),
                'n': len(data[cls]),
            }
        print(f'\n  {score_col} ({cell_grp}): KW H={H:.2f}, p={p_kw:.3e}')
        for pair, pval in pairs.items():
            print(f'    {pair}: p={pval:.3e}')

# ─────────────────────────────────────────────────────────────────
# 4. Myeloid subtype fractions by Class × Patient
# ─────────────────────────────────────────────────────────────────
print('\n[4/5] Myeloid fraction analysis ...')
# SPP1+ subtypes (key Pre-driver candidates)
SPP1_SUBTYPES = ['SPP1+A', 'SPP1+B']
SUPP_SUBTYPES = ['Anti-inflammatory', 'SPP1+A', 'SPP1+B']

frac_rows = []
for (patient, cls), grp in my.obs.groupby(['Patient', 'Class']):
    n_total = len(grp)
    if n_total < 5:
        continue
    n_spp1 = (grp['Cell_subtype'].isin(SPP1_SUBTYPES)).sum()
    n_supp = (grp['Cell_subtype'].isin(SUPP_SUBTYPES)).sum()
    n_prof = (grp['Cell_subtype'] == 'Pro-inflammatory').sum()
    frac_rows.append({
        'Patient': patient, 'Class': cls,
        'n_myeloid': n_total,
        'frac_SPP1': n_spp1 / n_total,
        'frac_suppressive': n_supp / n_total,
        'frac_proinflam': n_prof / n_total,
        'mean_Predriver_v2': grp['Predriver_v2'].mean(),
        'mean_TAM_score': grp['TAM_score'].mean(),
    })

frac_df = pd.DataFrame(frac_rows)
frac_df.to_csv(f'{RES}/myeloid_fractions_by_patient.csv', index=False)
print(frac_df.groupby('Class')[['frac_SPP1', 'frac_suppressive', 'mean_Predriver_v2']].mean().to_string())

# Stats: SPP1 fraction by Class
spp1_by_class = {cls: frac_df.loc[frac_df['Class'] == cls, 'frac_SPP1'].values
                 for cls in CLASS_ORDER}
H, p_kw, pairs = kruskal_pairwise(spp1_by_class)
print(f'\n  SPP1+ fraction by Class: KW H={H:.2f}, p={p_kw:.3e}')
for pair, pval in pairs.items():
    print(f'    {pair}: p={pval:.3e}')

# Infiltration index: ratio of suppressive to total myeloid per patient
infiltration = frac_df[['Patient', 'Class', 'frac_suppressive', 'frac_SPP1',
                         'mean_Predriver_v2', 'n_myeloid']].copy()
infiltration.to_csv(f'{RES}/infiltration_scores_v2.csv', index=False)

# ─────────────────────────────────────────────────────────────────
# 5. Plots
# ─────────────────────────────────────────────────────────────────
print('\n[5/5] Plotting ...')

# ── Fig A: Score distributions by Class (all cells) ───────────
score_cols = ['TAM_score', 'Predriver_v2', 'CD8Tex_score', 'Inflam_score']
score_labels = ['TAM score', 'Pre-driver score', 'CD8 Tex score', 'Inflammatory score']

fig, axes = plt.subplots(2, 2, figsize=(12, 9))
fig.suptitle('Immune cell scores by tumor region (GSE144735)',
             fontsize=12, fontweight='bold', fontfamily='Arial')
axes = axes.ravel()

for ax, score, label in zip(axes, score_cols, score_labels):
    data_by_class = [adata.obs.loc[adata.obs['Class'] == cls, score].values
                     for cls in CLASS_ORDER]
    vp = ax.violinplot(data_by_class, positions=range(3), showmedians=True)
    for pc, cls in zip(vp['bodies'], CLASS_ORDER):
        pc.set_facecolor(CLASS_COLORS[cls])
        pc.set_alpha(0.75)
    vp['cmedians'].set_color('black')
    vp['cmedians'].set_linewidth(2)
    ax.set_xticks(range(3))
    ax.set_xticklabels(CLASS_ORDER, fontsize=10)
    ax.set_ylabel(label, fontsize=9)
    ax.set_title(label, fontsize=10, fontweight='bold')

fig.tight_layout()
for ext in ('png', 'pdf'):
    path = f'{FIG}/score_by_class.{ext}'
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f'  Saved: {path}')
plt.close()

# ── Fig B: Myeloid subtype composition by Class ────────────────
subtype_order = ['Pro-inflammatory', 'Anti-inflammatory', 'SPP1+A', 'SPP1+B',
                 'Proliferating', 'cDC', 'Unknown']
subtype_colors = {
    'Pro-inflammatory':   '#4DBBD5',
    'Anti-inflammatory':  '#00A087',
    'SPP1+A':             '#F39B7F',
    'SPP1+B':             '#E64B35',
    'Proliferating':      '#8491B4',
    'cDC':                '#91D1C2',
    'Unknown':            '#B8B8B8',
}

fig, ax = plt.subplots(figsize=(8, 5))
my_counts = pd.crosstab(my.obs['Class'], my.obs['Cell_subtype']).reindex(CLASS_ORDER)
my_frac   = my_counts.div(my_counts.sum(axis=1), axis=0)

# Reindex subtype columns
plot_cols = [c for c in subtype_order if c in my_frac.columns]
my_frac_plot = my_frac[plot_cols]

bottom = np.zeros(len(CLASS_ORDER))
for st in plot_cols:
    vals = my_frac_plot[st].values
    ax.bar(CLASS_ORDER, vals, bottom=bottom,
           color=subtype_colors.get(st, '#CCCCCC'),
           label=st, edgecolor='white', linewidth=0.5)
    bottom += vals

ax.set_ylabel('Fraction', fontsize=10)
ax.set_title('Myeloid subtype composition by region', fontsize=11, fontweight='bold')
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)
ax.set_ylim(0, 1)
fig.tight_layout()
for ext in ('png', 'pdf'):
    path = f'{FIG}/myeloid_subtype_composition.{ext}'
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f'  Saved: {path}')
plt.close()

# ── Fig C: SPP1+ fraction & Pre-driver score per patient ──────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('Pre-driver myeloid enrichment at tumor border (GSE144735)',
             fontsize=11, fontweight='bold', fontfamily='Arial')

for ax, ycol, ylabel in [
    (axes[0], 'frac_SPP1', 'SPP1+ myeloid fraction'),
    (axes[1], 'mean_Predriver_v2', 'Mean Pre-driver score'),
]:
    for i, cls in enumerate(CLASS_ORDER):
        sub = frac_df[frac_df['Class'] == cls]
        xjit = np.random.default_rng(i).uniform(-0.2, 0.2, len(sub))
        ax.scatter(np.full(len(sub), i) + xjit, sub[ycol],
                   c=CLASS_COLORS[cls], s=60, alpha=0.8, zorder=3)
        # Boxplot overlay
        if len(sub) >= 3:
            bx = ax.boxplot(sub[ycol].values, positions=[i], widths=0.35,
                            patch_artist=True, showfliers=False)
            bx['boxes'][0].set_facecolor(CLASS_COLORS[cls])
            bx['boxes'][0].set_alpha(0.3)
    ax.set_xticks(range(3))
    ax.set_xticklabels(CLASS_ORDER, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(ylabel, fontsize=10, fontweight='bold')

fig.tight_layout()
for ext in ('png', 'pdf'):
    path = f'{FIG}/predriver_enrichment.{ext}'
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f'  Saved: {path}')
plt.close()

# ── Fig D: Cell type fractions by Class ───────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
ct_order = ['B cells', 'T cells', 'Myeloids', 'Stromal cells', 'Epithelial cells', 'Mast cells']
ct_colors = {'B cells': '#4DBBD5', 'T cells': '#E64B35', 'Myeloids': '#F39B7F',
             'Stromal cells': '#00A087', 'Epithelial cells': '#8491B4', 'Mast cells': '#91D1C2'}

all_ct = pd.crosstab(adata.obs['Class'], adata.obs['Cell_type']).reindex(CLASS_ORDER)
all_frac = all_ct.div(all_ct.sum(axis=1), axis=0)
plot_ct = [c for c in ct_order if c in all_frac.columns]

bottom = np.zeros(len(CLASS_ORDER))
for ct in plot_ct:
    vals = all_frac[ct].values
    ax.bar(CLASS_ORDER, vals, bottom=bottom,
           color=ct_colors.get(ct, '#CCCCCC'),
           label=ct, edgecolor='white', linewidth=0.5)
    bottom += vals

ax.set_ylabel('Fraction', fontsize=10)
ax.set_title('Cell type composition by tumor region', fontsize=11, fontweight='bold')
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)
ax.set_ylim(0, 1)
fig.tight_layout()
for ext in ('png', 'pdf'):
    path = f'{FIG}/celltype_composition.{ext}'
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f'  Saved: {path}')
plt.close()

# ── Fig E: Pre-driver score heatmap by Patient × Class ────────
pivot = frac_df.pivot_table(values='mean_Predriver_v2', index='Patient', columns='Class')
pivot = pivot.reindex(columns=CLASS_ORDER)
pivot = pivot.dropna(how='all')

fig, ax = plt.subplots(figsize=(6, max(4, len(pivot) * 0.4)))
im = ax.imshow(pivot.values, aspect='auto', cmap='YlOrRd',
               vmin=pivot.values[~np.isnan(pivot.values)].min(),
               vmax=pivot.values[~np.isnan(pivot.values)].max())
ax.set_xticks(range(len(CLASS_ORDER)))
ax.set_xticklabels(CLASS_ORDER, fontsize=10)
ax.set_yticks(range(len(pivot)))
ax.set_yticklabels(pivot.index, fontsize=8)
plt.colorbar(im, ax=ax, shrink=0.8, label='Mean Pre-driver score')
ax.set_title('Pre-driver score by Patient × Region', fontsize=10, fontweight='bold')
fig.tight_layout()
for ext in ('png', 'pdf'):
    path = f'{FIG}/predriver_heatmap_patient.{ext}'
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f'  Saved: {path}')
plt.close()

print('\n' + '=' * 65)
print('Summary — Border vs Core (GSE144735)')
print('=' * 65)
print(f'  Total cells: {adata.n_obs:,}')
print(f'  Myeloid cells: {my.n_obs:,}')
print(f'  Classes: {CLASS_ORDER}')
print(f'\n  SPP1+ Myeloid counts:')
print(my.obs[my.obs['Cell_subtype'].isin(SPP1_SUBTYPES)].groupby('Class').size().to_string())
print(f'\n  Results in : {RES}/')
print(f'  Figures in : {FIG}/')
