#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 1 Panel F — UMAP comparison: Unintegrated (PCA) / Harmony / scVI
2 rows × 3 columns:
  Row 1: coloured by dataset (batch)
  Row 2: coloured by cell type
"""
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os, sys, io
from datetime import datetime

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

os.makedirs('figures/phase1', exist_ok=True)

print("=" * 60)
print(f"Figure 1F — UMAP Comparison  |  {datetime.now().strftime('%H:%M:%S')}")
print("=" * 60)

# ─────────────────────────────────────────────────────────────────
# 1. Load data
# ─────────────────────────────────────────────────────────────────
print("\n[1/4] Loading data ...")
adata = sc.read_h5ad('processed_data/merged_annotated.h5ad')
print(f"  {adata.n_obs:,} cells  |  X_scvi: {adata.obsm['X_scvi'].shape}")

batch_labels    = adata.obs['batch'].astype(str).values
celltype_labels = adata.obs['cell_type_grouped'].astype(str).values

BATCH_COLORS = {
    'GSE132257': '#E64B35',
    'GSE132465': '#4DBBD5',
    'GSE125527': '#00A087',
    'GSE150115': '#3C5488',
    'cellxgene': '#F39B7F',
}
TYPE_COLORS = {
    'T cells'   : '#E64B35',
    'NK cells'  : '#3C5488',
    'B cells'   : '#4DBBD5',
    'Myeloids'  : '#00A087',
    'Mast cells': '#F39B7F',
    'Epithelial': '#91D1C2',
    'Stromal'   : '#DC0000',
    'Malignant' : '#B09C85',
    'Unknown'   : '#AAAAAA',
}

# ─────────────────────────────────────────────────────────────────
# 2. PCA embedding → UMAP
# ─────────────────────────────────────────────────────────────────
print("\n[2/4] PCA + UMAP (unintegrated) ...")
adata_pca = adata.copy()
sc.pp.normalize_total(adata_pca, target_sum=1e4)
sc.pp.log1p(adata_pca)
sc.pp.highly_variable_genes(adata_pca, n_top_genes=3000, batch_key='batch')
sc.pp.scale(adata_pca, max_value=10)
sc.tl.pca(adata_pca, n_comps=30, use_highly_variable=True)
sc.pp.neighbors(adata_pca, use_rep='X_pca', n_neighbors=30)
sc.tl.umap(adata_pca, min_dist=0.3)
umap_pca = adata_pca.obsm['X_umap']
X_pca    = adata_pca.obsm['X_pca']
del adata_pca
print(f"  PCA UMAP shape: {umap_pca.shape}")

# ─────────────────────────────────────────────────────────────────
# 3. Harmony embedding → UMAP
# ─────────────────────────────────────────────────────────────────
print("\n[3/4] Harmony + UMAP ...")
import harmonypy as hm
ho = hm.run_harmony(
    X_pca.astype(np.float64),
    pd.DataFrame({'batch': batch_labels}),
    'batch',
    max_iter_harmony=20,
    random_state=0,
    verbose=True,
)
X_harmony = np.array(ho.Z_corr)

# Build a lightweight AnnData just for neighbors/UMAP
adata_hm = sc.AnnData(obs=adata.obs.copy())
adata_hm.obsm['X_harmony'] = X_harmony.astype(np.float32)
sc.pp.neighbors(adata_hm, use_rep='X_harmony', n_neighbors=30)
sc.tl.umap(adata_hm, min_dist=0.3)
umap_hm = adata_hm.obsm['X_umap']
del adata_hm
print(f"  Harmony UMAP shape: {umap_hm.shape}")

# ─────────────────────────────────────────────────────────────────
# 4. scVI UMAP (pre-computed in merged_annotated.h5ad)
# ─────────────────────────────────────────────────────────────────
umap_scvi = adata.obsm['X_umap']
print(f"\n  scVI UMAP shape: {umap_scvi.shape}  (pre-computed)")

# ─────────────────────────────────────────────────────────────────
# 5. Plot 2 × 3 grid
# ─────────────────────────────────────────────────────────────────
print("\n[4/4] Plotting ...")

METHODS = [
    ('Unintegrated (PCA)', umap_pca),
    ('Harmony',            umap_hm),
    ('scVI',               umap_scvi),
]

fig, axes = plt.subplots(2, 3, figsize=(14, 9))
fig.suptitle('Figure 1F — UMAP Comparison: Unintegrated / Harmony / scVI',
             fontsize=13, fontweight='bold', y=1.01)

DOT_SIZE  = 0.4
DOT_ALPHA = 0.35

def scatter_umap(ax, xy, labels, color_map, title, subtitle=None):
    """Plot UMAP scatter coloured by category."""
    # Draw in consistent order so no group hides another; shuffle first
    rng = np.random.default_rng(0)
    idx = rng.permutation(len(labels))
    xy_s, lab_s = xy[idx], labels[idx]

    for label, color in color_map.items():
        mask = lab_s == label
        if mask.sum() == 0:
            continue
        ax.scatter(xy_s[mask, 0], xy_s[mask, 1],
                   c=color, s=DOT_SIZE, alpha=DOT_ALPHA,
                   linewidths=0, rasterized=True)
    if subtitle:
        ax.set_title(f'{title}\n{subtitle}', fontsize=10, fontweight='bold', pad=4)
    else:
        ax.set_title(title, fontsize=10, fontweight='bold', pad=4)
    ax.axis('off')

# ── Row 0: coloured by batch ─────────────────────────────────────
for col, (method, umap_xy) in enumerate(METHODS):
    scatter_umap(axes[0, col], umap_xy, batch_labels, BATCH_COLORS,
                 method, subtitle='by dataset')

# ── Row 1: coloured by cell type ─────────────────────────────────
for col, (method, umap_xy) in enumerate(METHODS):
    scatter_umap(axes[1, col], umap_xy, celltype_labels, TYPE_COLORS,
                 method, subtitle='by cell type')

# ── Legends (right margin of each row) ───────────────────────────
# Batch legend: attach to axes[0, 2]
batch_handles = [mpatches.Patch(facecolor=c, label=b, linewidth=0)
                 for b, c in BATCH_COLORS.items()]
axes[0, 2].legend(handles=batch_handles, title='Dataset',
                  title_fontsize=8, fontsize=7.5,
                  loc='lower right', framealpha=0.85,
                  edgecolor='#cccccc', markerscale=1.5,
                  handlelength=1.2, borderpad=0.6)

# Cell type legend: attach to axes[1, 2]
type_handles = [mpatches.Patch(facecolor=c, label=t, linewidth=0)
                for t, c in TYPE_COLORS.items()
                if t in np.unique(celltype_labels)]
axes[1, 2].legend(handles=type_handles, title='Cell type',
                  title_fontsize=8, fontsize=7.5,
                  loc='lower right', framealpha=0.85,
                  edgecolor='#cccccc', markerscale=1.5,
                  handlelength=1.2, borderpad=0.6)

plt.tight_layout()

out_png = 'figures/phase1/figure1_panel_F_umap.png'
out_pdf = 'figures/phase1/figure1_panel_F_umap.pdf'
plt.savefig(out_png, dpi=200, bbox_inches='tight', facecolor='white')
plt.savefig(out_pdf, bbox_inches='tight', facecolor='white')
plt.close()
print(f"\nSaved: {out_png}")
print(f"Saved: {out_pdf}")
print("Done.")
