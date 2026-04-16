#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CellAssign cell type annotation for merged scVI dataset.
Uses 40 marker genes confirmed present in the 3k HVG set.
Annotates all 195,090 cells across all 5 batches.
"""

import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scvi
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os, warnings, sys, io
from datetime import datetime

warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
plt.rcParams['font.family'] = 'Arial'
os.makedirs('figures/phase1', exist_ok=True)

# ── Load ──────────────────────────────────────────────────────────────────────
print(f"[{datetime.now():%H:%M:%S}] Loading data...")
adata = sc.read_h5ad('scvi_results_hvg3k/merged_scvi_hvg3k.h5ad')
print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# ── Marker gene matrix ────────────────────────────────────────────────────────
# Binary: 1 = marker for this cell type
# All genes below confirmed present in HVG3k set

CELL_TYPES = [
    'CD4 T', 'CD8 T', 'Treg', 'NK',
    'B cell', 'Plasma',
    'Myeloid', 'Mast',
    'Stromal', 'Epithelial',
]

MARKERS = {
    'CD4 T':     ['CD3D', 'CD3E', 'CD4', 'IL7R', 'CCR7'],
    'CD8 T':     ['CD3D', 'CD3E', 'CD8A', 'GZMA', 'GZMB', 'GZMK', 'NKG7'],
    'Treg':      ['CD3D', 'CD3E', 'FOXP3', 'IL2RA', 'IKZF2', 'CTLA4'],
    'NK':        ['NCAM1', 'GNLY', 'NKG7', 'KLRB1', 'KLRC1', 'KLRD1', 'GZMB', 'PRF1'],
    'B cell':    ['CD19', 'CD79A', 'CD79B', 'MS4A1'],
    'Plasma':    ['MZB1', 'SDC1', 'IGHA1', 'IGHG1', 'IGKC', 'DERL3'],
    'Myeloid':   ['LYZ', 'CD14', 'CD68', 'CSF1R', 'FCGR3A', 'ITGAM'],
    'Mast':      ['FCER1A'],
    'Stromal':   ['VIM', 'ACTA2'],
    'Epithelial':['KRT18'],
}

# Verify all genes are available
avail = set(adata.var_names)
for ct, genes in MARKERS.items():
    missing = [g for g in genes if g not in avail]
    if missing:
        print(f"  WARNING: {ct} missing: {missing}")
        MARKERS[ct] = [g for g in genes if g in avail]

all_markers = sorted(set(g for gs in MARKERS.values() for g in gs))
print(f"\nUsing {len(all_markers)} marker genes: {all_markers}")

# Build binary DataFrame: index=genes, columns=cell_types
marker_df = pd.DataFrame(0, index=all_markers, columns=CELL_TYPES)
for ct, genes in MARKERS.items():
    for g in genes:
        if g in marker_df.index:
            marker_df.loc[g, ct] = 1

print("\nMarker gene matrix (transposed for readability):")
print(marker_df.T.to_string())

# ── Compute size factor from full counts (before subsetting) ──────────────────
print(f"\n[{datetime.now():%H:%M:%S}] Computing size factors...")
if 'counts' in adata.layers:
    counts_full = adata.layers['counts']
else:
    counts_full = adata.X
if sp.issparse(counts_full):
    lib_size = np.asarray(counts_full.sum(axis=1)).flatten()
else:
    lib_size = counts_full.sum(axis=1)
# Normalize: divide by median library size (standard CellAssign convention)
size_factor = lib_size / np.median(lib_size)
print(f"  Size factor range: [{size_factor.min():.3f}, {size_factor.max():.3f}], median=1.000")

# ── Subset adata to marker genes ──────────────────────────────────────────────
print(f"[{datetime.now():%H:%M:%S}] Subsetting to {len(all_markers)} marker genes...")
adata_ca = adata[:, all_markers].copy()

# CellAssign requires raw counts in X
if 'counts' in adata.layers:
    adata_ca.X = adata[:, all_markers].layers['counts'].copy()
# (already raw if using counts layer)

adata_ca.obs['size_factor'] = size_factor

# ── Setup & train CellAssign ──────────────────────────────────────────────────
print(f"[{datetime.now():%H:%M:%S}] Setting up CellAssign...")
scvi.external.CellAssign.setup_anndata(adata_ca, size_factor_key='size_factor')

model = scvi.external.CellAssign(adata_ca, cell_type_markers=marker_df)
print(f"[{datetime.now():%H:%M:%S}] Training CellAssign (max_epochs=200)...")
model.train(max_epochs=200)

# ── Predict ───────────────────────────────────────────────────────────────────
print(f"[{datetime.now():%H:%M:%S}] Predicting cell types...")
predictions = model.predict()   # DataFrame: cells × cell_types (probabilities)

adata.obs['cellass_type'] = predictions.idxmax(axis=1).values
adata.obs['cellass_prob_max'] = predictions.max(axis=1).values
for ct in CELL_TYPES:
    adata.obs[f'cellass_prob_{ct}'] = predictions[ct].values

# ── Summary ───────────────────────────────────────────────────────────────────
print("\n=== CellAssign cell type distribution ===")
ct_counts = adata.obs['cellass_type'].value_counts()
for ct, n in ct_counts.items():
    pct = n / adata.n_obs * 100
    print(f"  {ct:15s}: {n:7,}  ({pct:5.1f}%)")

print("\n=== Per-batch breakdown (%) ===")
batch_ct = pd.crosstab(adata.obs['batch'], adata.obs['cellass_type'],
                       normalize='index') * 100
ct_order_print = [c for c in CELL_TYPES if c in batch_ct.columns]
print(batch_ct[ct_order_print].round(1).to_string())

print(f"\nMedian confidence (max prob): {adata.obs['cellass_prob_max'].median():.3f}")

# ── Save updated h5ad ─────────────────────────────────────────────────────────
out_path = 'scvi_results_hvg3k/merged_scvi_hvg3k.h5ad'
print(f"\n[{datetime.now():%H:%M:%S}] Saving {out_path}...")
adata.write_h5ad(out_path)
print("  Saved.")

# ── UMAP plots ────────────────────────────────────────────────────────────────
CT_COLORS = {
    'CD4 T':     '#4DBBD5',
    'CD8 T':     '#00A087',
    'Treg':      '#3C5488',
    'NK':        '#F39B7F',
    'B cell':    '#F7DC6F',
    'Plasma':    '#8491B4',
    'Myeloid':   '#91D1C2',
    'Mast':      '#E64B35',
    'Stromal':   '#B09C85',
    'Epithelial':'#BBCC33',
}

xy  = adata.obsm['X_umap']
rng = np.random.default_rng(42)
idx = rng.permutation(adata.n_obs)
ct_vals = adata.obs['cellass_type'].values[idx]

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
fig.subplots_adjust(left=0.04, right=0.98, top=0.92, bottom=0.08, wspace=0.30)

# ── Left: UMAP by CellAssign type ─────────────────────────────────────────────
ax = axes[0]
for ct in CELL_TYPES:
    m = ct_vals == ct
    if m.sum() == 0:
        continue
    ax.scatter(xy[idx[m], 0], xy[idx[m], 1],
               c=CT_COLORS.get(ct, '#aaa'), s=0.3, alpha=0.45,
               linewidths=0, rasterized=True)

leg = [Patch(color=CT_COLORS[ct], label=ct) for ct in CELL_TYPES if ct in ct_counts.index]
ax.legend(handles=leg, fontsize=7.5, frameon=False, ncol=1,
          loc='lower right', handlelength=0.9, handleheight=0.9,
          borderpad=0.3, labelspacing=0.35)
ax.set_xlabel('UMAP 1', fontsize=9)
ax.set_ylabel('UMAP 2', fontsize=9)
ax.tick_params(labelsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_title('CellAssign Annotation', fontsize=10, pad=4)
ax.text(0.02, 0.98, f'n={adata.n_obs:,}', transform=ax.transAxes,
        fontsize=7, va='top', color='#555')

# ── Right: per-batch stacked bar ───────────────────────────────────────────────
ax = axes[1]
batches_plot = ['GSE132257', 'GSE132465', 'cellxgene', 'GSE125527', 'GSE150115']
x  = np.arange(len(batches_plot))
W  = 0.72
bottoms = np.zeros(len(batches_plot))

for ct in CELL_TYPES:
    vals = np.array([
        batch_ct.loc[b, ct] if b in batch_ct.index and ct in batch_ct.columns else 0.0
        for b in batches_plot
    ])
    ax.bar(x, vals, bottom=bottoms, color=CT_COLORS[ct],
           width=W, linewidth=0.3, edgecolor='white')
    bottoms += vals

# Cell counts below
for xi, b in enumerate(batches_plot):
    n = (adata.obs['batch'] == b).sum()
    ax.text(xi, -4, f'n={n:,}', ha='center', va='top', fontsize=6.5, color='#444')

ax.set_xticks(x)
ax.set_xticklabels(['GSE132257', 'GSE132465', 'CXG', 'GSE125527', 'GSE150115'],
                   fontsize=8)
ax.set_ylabel('Proportion (%)', fontsize=9)
ax.set_xlim(-0.6, 4.6)
ax.set_ylim(-10, 108)
ax.tick_params(labelsize=7.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_title('Cell Type Composition per Batch (CellAssign)', fontsize=10, pad=4)

leg2 = [Patch(color=CT_COLORS[ct], label=ct) for ct in CELL_TYPES]
ax.legend(handles=leg2, fontsize=7.5, frameon=False, ncol=1,
          loc='center left', bbox_to_anchor=(1.01, 0.5),
          handlelength=0.9, handleheight=0.9, labelspacing=0.4)

for ext in ['png', 'pdf']:
    fig.savefig(f'figures/phase1/cellass_annotation.{ext}',
                dpi=200, bbox_inches='tight')
plt.close(fig)
print(f"[{datetime.now():%H:%M:%S}] Saved figures/phase1/cellass_annotation.png/.pdf")
print("Done.")
