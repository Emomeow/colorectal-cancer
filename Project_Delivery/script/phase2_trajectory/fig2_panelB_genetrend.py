#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2 – Panel B: Gene expression heatmap along pseudotime
Uses top pseudotime-correlated genes; loads expression from h5ad
in backed mode, bins cells by pseudotime, smooths, and plots heatmap.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scanpy as sc
import os
from scipy.ndimage import uniform_filter1d

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
plt.rcParams['font.family'] = 'Arial'

RES = 'phase2_results/myeloid_v3'
FIG_OUT = 'figures/phase2/myeloid_v3'
os.makedirs(FIG_OUT, exist_ok=True)

# ── Load pseudotime + gene correlation ────────────────────────────
meta = pd.read_csv(f'{RES}/cell_metadata.csv', index_col=0)
pt_corr = pd.read_csv(f'{RES}/pseudotime_gene_correlation.csv')

# Select top 8 rising (positive) and 8 falling (negative) genes
rising  = pt_corr[pt_corr['rho'] > 0].head(8)['gene'].tolist()
falling = pt_corr[pt_corr['rho'] < 0].head(8)['gene'].tolist()
genes_sel = falling + rising
print(f'Falling: {falling}')
print(f'Rising:  {rising}')

# ── Load h5ad subset (backed) for expression ──────────────────────
print('Loading h5ad (backed, row-only subset)...')
adata_all = sc.read_h5ad('processed_data/merged_annotated_fine.h5ad', backed='r')
mask = (adata_all.obs['batch'].isin(['GSE125527', 'GSE132465'])) & \
       (adata_all.obs['cell_type_grouped'] == 'Myeloids')
# Step 1: load only the myeloid rows into memory (row indexing only)
adata_sub = adata_all[mask].to_memory()
adata_all.file.close()
del adata_all
# Step 2: then subset to selected genes in memory
avail_genes = [g for g in genes_sel if g in adata_sub.var_names]
adata_sub = adata_sub[:, avail_genes].copy()

print(f'Loaded {adata_sub.n_obs} cells x {adata_sub.n_vars} genes')
genes_sel = [g for g in genes_sel if g in adata_sub.var_names]
falling   = [g for g in falling   if g in genes_sel]
rising    = [g for g in rising    if g in genes_sel]
genes_ordered = falling + rising

# ── Attach pseudotime ─────────────────────────────────────────────
common_idx = adata_sub.obs_names.intersection(meta.index)
adata_sub  = adata_sub[common_idx]
pt = meta.loc[common_idx, 'palantir_pseudotime'].values
adata_sub.obs['pseudotime'] = pt

# ── Sort by pseudotime, bin into 100 windows, take mean ──────────
N_BINS = 100
order  = np.argsort(pt)
chunk  = max(1, len(order) // N_BINS)

X = adata_sub[:, genes_ordered].X
if hasattr(X, 'toarray'):
    X = X.toarray()

bin_expr = []
bin_pt   = []
for i in range(0, len(order), chunk):
    idx = order[i:i+chunk]
    bin_expr.append(X[idx].mean(axis=0))
    bin_pt.append(pt[idx].mean())

mat = np.array(bin_expr).T   # genes × bins
pt_bins = np.array(bin_pt)

# Smooth each gene row
SMOOTH = 8
mat_sm = uniform_filter1d(mat, size=SMOOTH, axis=1)

# Z-score per gene
mat_z = (mat_sm - mat_sm.mean(axis=1, keepdims=True)) / \
        (mat_sm.std(axis=1, keepdims=True) + 1e-9)
mat_z = np.clip(mat_z, -2.5, 2.5)

# ── Pre-driver window (branch point-based, data-driven) ───────────
# Load branch probabilities to compute branch point
bp_df = pd.read_csv(f'{RES}/branch_probs.csv', index_col=0)
common_bp = bp_df.index.intersection(meta.index)
bp_vals = bp_df.loc[common_bp].values
pt_for_bp = meta.loc[common_bp, 'palantir_pseudotime'].values

# Branch point = pseudotime at max fate-prob variance
n_bins_bp = 50
bins_bp   = np.linspace(pt_for_bp.min(), pt_for_bp.max(), n_bins_bp + 1)
bin_idx   = np.clip(np.digitize(pt_for_bp, bins_bp) - 1, 0, n_bins_bp - 1)
variances = [bp_vals[bin_idx == b].var(axis=0).mean()
             if (bin_idx == b).sum() > 2 else 0
             for b in range(n_bins_bp)]
branch_pt = ((bins_bp[:-1] + bins_bp[1:]) / 2)[np.argmax(variances)]

pt_span = pt_for_bp.max() - pt_for_bp.min()
PD_LO   = branch_pt - 0.30 * pt_span   # 30% span before branch point
PD_HI   = branch_pt                    # upper boundary = branch point
print(f'Branch point: {branch_pt:.3f}  Pre-driver window: [{PD_LO:.3f}, {PD_HI:.3f}]')

pd_start = np.searchsorted(pt_bins, PD_LO)
pd_end   = np.searchsorted(pt_bins, PD_HI)

# ── Plot ──────────────────────────────────────────────────────────
n_genes   = len(genes_ordered)
n_falling = len(falling)
n_rising  = len(rising)

fig, ax = plt.subplots(figsize=(11, max(4.5, n_genes * 0.42 + 1.5)))
fig.subplots_adjust(left=0.18, right=0.88, top=0.89, bottom=0.12)

im = ax.imshow(mat_z, aspect='auto', cmap='RdBu_r',
               vmin=-2.5, vmax=2.5, interpolation='nearest')

# Pre-driver window
ax.axvspan(pd_start - 0.5, pd_end + 0.5, color='#FFD700', alpha=0.22)
ax.axvline(pd_start - 0.5, color='#B8860B', lw=1.2, ls='--')
ax.axvline(pd_end   + 0.5, color='#B8860B', lw=1.2, ls='--')

# Group separator
if n_falling > 0 and n_rising > 0:
    ax.axhline(n_falling - 0.5, color='white', lw=1.5)

# Gene labels with colors
ax.set_yticks(range(n_genes))
gene_cols = ['#4DBBD5'] * n_falling + ['#E64B35'] * n_rising
ax.set_yticklabels(genes_ordered, fontsize=9)
for lbl, col in zip(ax.get_yticklabels(), gene_cols):
    lbl.set_color(col)

# Pseudotime x-axis
tick_vals = [0, 0.25, 0.5, 0.75, 1.0]
tick_pos  = [np.searchsorted(pt_bins, v) for v in tick_vals]
ax.set_xticks(tick_pos)
ax.set_xticklabels([f'{v:.2f}' for v in tick_vals], fontsize=9)
ax.set_xlabel('Pseudotime', fontsize=10)

ax.set_title(
    'Key gene dynamics along myeloid trajectory\n'
    '(Classical_Mono → Suppressive TAM, CRC lineage)',
    fontsize=10, pad=6
)

cbar = plt.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
cbar.set_label('Z-score', fontsize=8.5)
cbar.ax.tick_params(labelsize=8)

rising_p  = mpatches.Patch(color='#E64B35', label='Rising (suppressive)')
falling_p = mpatches.Patch(color='#4DBBD5', label='Falling (inflammatory)')
pd_p      = mpatches.Patch(color='#FFD700', alpha=0.7, label='Pre-driver window')
ax.legend(handles=[falling_p, rising_p, pd_p],
          fontsize=8, loc='lower right', frameon=True, framealpha=0.85)

ax.text(-0.14, 1.04, 'B', transform=ax.transAxes,
        fontsize=14, fontweight='bold', va='bottom')

for ext in ('png', 'pdf'):
    path = f'{FIG_OUT}/gene_trend_heatmap_panelB.{ext}'
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f'Saved: {path}')
plt.close()
print('Done.')
