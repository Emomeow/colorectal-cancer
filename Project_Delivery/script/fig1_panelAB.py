#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 1 – Panel A, B & C: UMAP colored by batch / disease / cell type"""

import scanpy as sc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
plt.rcParams['font.family'] = 'Arial'
os.makedirs('figures/phase1', exist_ok=True)

adata = sc.read_h5ad('scvi_results_hvg3k/merged_scvi_hvg3k.h5ad')
xy = adata.obsm['X_umap']

BATCH_COLORS = {
    'GSE132257': '#E64B35',
    'GSE132465': '#4DBBD5',
    'GSE125527': '#00A087',
    'GSE150115': '#3C5488',
    'cellxgene': '#F39B7F',
}
DISEASE_COLORS = {'CRC': '#E64B35', 'UC': '#4DBBD5'}

BATCH_LABEL = {
    'GSE132257': 'GSE132257 (CRC)',
    'GSE132465': 'GSE132465 (CRC)',
    'GSE125527': 'GSE125527 (UC)',
    'GSE150115': 'GSE150115 (UC)',
    'cellxgene': 'cellxgene (CRC)',
}

CT_ORDER = ['CD4 T', 'CD8 T', 'Treg', 'NK',
            'B cell', 'Plasma', 'Myeloid', 'Mast', 'Stromal', 'Epithelial']
CT_COLORS = {
    'CD4 T':     '#E64B35', 'CD8 T':  '#4DBBD5', 'Treg':    '#00A087',
    'NK':        '#3C5488', 'B cell': '#F39B7F',  'Plasma':  '#8491B4',
    'Myeloid':   '#91D1C2', 'Mast':   '#DC0000',  'Stromal': '#7E6148',
    'Epithelial':'#B09C85',
}

fig, axes = plt.subplots(1, 3, figsize=(17, 5.5))
fig.subplots_adjust(left=0.03, right=0.98, top=0.92, bottom=0.06, wspace=0.28)

# ── Panel A: batch ────────────────────────────────────────────────────────────
ax = axes[0]
batches = sorted(adata.obs['batch'].unique())
# Shuffle so no one batch is always on top
rng = np.random.default_rng(42)
idx = rng.permutation(adata.n_obs)
batch_vals = adata.obs['batch'].values[idx]
for b in batches:
    m = batch_vals == b
    ax.scatter(xy[idx[m], 0], xy[idx[m], 1],
               c=BATCH_COLORS[b], s=0.3, alpha=0.4,
               linewidths=0, rasterized=True)

leg = [Patch(color=BATCH_COLORS[b], label=BATCH_LABEL[b]) for b in batches]
ax.legend(handles=leg, fontsize=7.5, frameon=False, ncol=1,
          loc='lower right', handlelength=0.9, handleheight=0.9,
          borderpad=0.3, labelspacing=0.35)
ax.set_xlabel('UMAP 1', fontsize=9)
ax.set_ylabel('UMAP 2', fontsize=9)
ax.tick_params(labelsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.text(-0.06, 1.02, 'A', transform=ax.transAxes,
        fontsize=14, fontweight='bold', va='bottom')
ax.set_title('Dataset', fontsize=10, pad=4)

# ── Panel B: disease ──────────────────────────────────────────────────────────
ax = axes[1]
dis_vals = adata.obs['disease'].values[idx]
for dis in ['UC', 'CRC']:          # UC first so CRC is on top (smaller set)
    m = dis_vals == dis
    ax.scatter(xy[idx[m], 0], xy[idx[m], 1],
               c=DISEASE_COLORS[dis], s=0.3, alpha=0.45,
               linewidths=0, rasterized=True)

leg = [Patch(color=DISEASE_COLORS[d], label=d) for d in ['CRC', 'UC']]
ax.legend(handles=leg, fontsize=8.5, frameon=False, loc='lower right',
          handlelength=0.9, handleheight=0.9)
ax.set_xlabel('UMAP 1', fontsize=9)
ax.set_ylabel('UMAP 2', fontsize=9)
ax.tick_params(labelsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.text(-0.06, 1.02, 'B', transform=ax.transAxes,
        fontsize=14, fontweight='bold', va='bottom')
ax.set_title('Disease', fontsize=10, pad=4)

# Cell counts annotation
for dis in ['CRC', 'UC']:
    n = (adata.obs['disease'] == dis).sum()
    ax.text(0.02, 0.98 if dis == 'CRC' else 0.91,
            f'{dis}: {n:,}', transform=ax.transAxes,
            fontsize=7, color=DISEASE_COLORS[dis], va='top')

# ── Panel C: cell type ────────────────────────────────────────────────────────
ax = axes[2]
ct_vals = adata.obs['cellass_type'].values[idx]
for ct in CT_ORDER:
    m = ct_vals == ct
    if m.sum() == 0:
        continue
    ax.scatter(xy[idx[m], 0], xy[idx[m], 1],
               c=CT_COLORS.get(ct, '#aaaaaa'), s=0.3, alpha=0.45,
               linewidths=0, rasterized=True)

leg = [Patch(color=CT_COLORS[ct], label=ct) for ct in CT_ORDER
       if ct in adata.obs['cellass_type'].values]
ax.legend(handles=leg, fontsize=6.5, frameon=False, ncol=1,
          loc='lower right', handlelength=0.9, handleheight=0.9,
          borderpad=0.3, labelspacing=0.3)
ax.set_xlabel('UMAP 1', fontsize=9)
ax.set_ylabel('UMAP 2', fontsize=9)
ax.tick_params(labelsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.text(-0.06, 1.02, 'C', transform=ax.transAxes,
        fontsize=14, fontweight='bold', va='bottom')
ax.set_title('Cell Type (CellAssign)', fontsize=10, pad=4)

for ext in ['png', 'pdf']:
    fig.savefig(f'figures/phase1/fig1_panelAB.{ext}',
                dpi=200, bbox_inches='tight')
plt.close(fig)
print("Saved figures/phase1/fig1_panelAB.png/.pdf")
