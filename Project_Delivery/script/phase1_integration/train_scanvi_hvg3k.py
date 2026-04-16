#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scANVI: Semi-supervised extension of scVI.
Loads pre-trained scVI model, fine-tunes with CellAssign cell type labels,
aiming to improve batch mixing (iLISI) while maintaining bio conservation.

Input:
  scvi_results_hvg3k/merged_scvi_hvg3k.h5ad  (has cellass_type + X_scvi)
  scvi_model_hvg3k/                           (pre-trained scVI model)

Output:
  scanvi_model_hvg3k/                         (scANVI model)
  scanvi_results_hvg3k/merged_scanvi_hvg3k.h5ad
  figures/phase1/fig_scanvi_umap.png/.pdf
"""

import scanpy as sc
import scvi
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os, sys, io, warnings
from datetime import datetime

warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
plt.rcParams['font.family'] = 'Arial'

DATA_PATH   = 'scvi_results_hvg3k/merged_scvi_hvg3k.h5ad'
SCVI_DIR    = 'scvi_model_hvg3k'
SCANVI_DIR  = 'scanvi_model_hvg3k'
RESULTS_DIR = 'scanvi_results_hvg3k'
FIG_DIR     = 'figures/phase1'

MAX_EPOCHS         = 50
EARLY_STOP_PATIENCE = 10
N_SAMPLES_PER_LABEL = 100   # labeled cells sampled per mini-batch

os.makedirs(SCANVI_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# ── [1] Load data ──────────────────────────────────────────────────────────────
print(f"[{datetime.now():%H:%M:%S}] Loading data...")
adata = sc.read_h5ad(DATA_PATH)
print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"  cellass_type: {adata.obs['cellass_type'].nunique()} types, "
      f"{(adata.obs['cellass_type'] == 'Unknown').sum()} unlabeled")

# ── [2] Setup scVI registry (must match training config exactly) ───────────────
print(f"\n[{datetime.now():%H:%M:%S}] Setting up scVI registry...")
scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    batch_key='batch',
    categorical_covariate_keys=['disease'],
)

# ── [3] Load pre-trained scVI model ───────────────────────────────────────────
print(f"\n[{datetime.now():%H:%M:%S}] Loading scVI model from {SCVI_DIR}/...")
scvi_model = scvi.model.SCVI.load(SCVI_DIR, adata=adata)
print(f"  Loaded: {scvi_model}")

# ── [4] Initialize scANVI from scVI ───────────────────────────────────────────
print(f"\n[{datetime.now():%H:%M:%S}] Initializing scANVI from scVI model...")
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    unlabeled_category='Unknown',   # cells with this label are treated as unlabeled
    labels_key='cellass_type',
)
print(f"  scANVI model: {scanvi_model.summary_stats}")

# ── [5] Train scANVI ──────────────────────────────────────────────────────────
print(f"\n[{datetime.now():%H:%M:%S}] Training scANVI "
      f"(max_epochs={MAX_EPOCHS}, patience={EARLY_STOP_PATIENCE})...")
start = datetime.now()
scanvi_model.train(
    max_epochs=MAX_EPOCHS,
    n_samples_per_label=N_SAMPLES_PER_LABEL,
    early_stopping=True,
    early_stopping_patience=EARLY_STOP_PATIENCE,
    early_stopping_monitor='elbo_validation',
)
duration = (datetime.now() - start).total_seconds()
print(f"  Finished in {duration/60:.1f} min")

# ── [6] Save scANVI model ─────────────────────────────────────────────────────
scanvi_model.save(SCANVI_DIR, overwrite=True)
print(f"  Saved → {SCANVI_DIR}/")

# ── [7] Get scANVI latent representation ──────────────────────────────────────
print(f"\n[{datetime.now():%H:%M:%S}] Getting latent representation...")
adata.obsm['X_scanvi'] = scanvi_model.get_latent_representation()
print(f"  X_scanvi shape: {adata.obsm['X_scanvi'].shape}")

# ── [8] Neighbors + UMAP on scANVI embedding ──────────────────────────────────
print(f"\n[{datetime.now():%H:%M:%S}] Building neighbors + UMAP (scANVI)...")

# Backup scVI UMAP
adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()

sc.pp.neighbors(adata, use_rep='X_scanvi', n_neighbors=15, key_added='scanvi_neighbors')
sc.tl.umap(adata, neighbors_key='scanvi_neighbors', min_dist=0.3)
adata.obsm['X_umap_scanvi'] = adata.obsm['X_umap'].copy()

# Restore scVI UMAP as default
adata.obsm['X_umap'] = adata.obsm['X_umap_scvi'].copy()
print(f"  X_umap_scanvi computed, X_umap restored to scVI embedding")

# ── [9] Save h5ad ─────────────────────────────────────────────────────────────
out_h5ad = f'{RESULTS_DIR}/merged_scanvi_hvg3k.h5ad'
adata.write_h5ad(out_h5ad, compression='gzip')
print(f"\n[{datetime.now():%H:%M:%S}] Saved: {out_h5ad} "
      f"({os.path.getsize(out_h5ad)/1024**2:.1f} MB)")

# ── [10] Comparison figure: scVI vs scANVI UMAP ───────────────────────────────
print(f"\n[{datetime.now():%H:%M:%S}] Plotting comparison figure...")

BATCH_COLORS = {
    'GSE132257': '#E64B35', 'GSE132465': '#4DBBD5',
    'GSE125527': '#00A087', 'GSE150115': '#3C5488', 'cellxgene': '#F39B7F',
}
DISEASE_COLORS = {'CRC': '#E64B35', 'UC': '#4DBBD5'}
CT_ORDER = ['CD4 T', 'CD8 T', 'Treg', 'NK',
            'B cell', 'Plasma', 'Myeloid', 'Mast', 'Stromal', 'Epithelial']
CT_COLORS = {
    'CD4 T': '#E64B35', 'CD8 T': '#4DBBD5', 'Treg': '#00A087',
    'NK': '#3C5488', 'B cell': '#F39B7F', 'Plasma': '#8491B4',
    'Myeloid': '#91D1C2', 'Mast': '#DC0000', 'Stromal': '#7E6148',
    'Epithelial': '#B09C85',
}

rng = np.random.default_rng(42)
idx = rng.permutation(adata.n_obs)
batch_vals = adata.obs['batch'].values[idx]
dis_vals   = adata.obs['disease'].values[idx]
ct_vals    = adata.obs['cellass_type'].values[idx]

fig, axes = plt.subplots(2, 3, figsize=(18, 11))
fig.subplots_adjust(left=0.04, right=0.97, top=0.93, bottom=0.06,
                    wspace=0.28, hspace=0.35)

embeddings = {
    'scVI':   adata.obsm['X_umap_scvi'][idx],
    'scANVI': adata.obsm['X_umap_scanvi'][idx],
}

for row, (method, xy) in enumerate(embeddings.items()):
    # Batch
    ax = axes[row, 0]
    for b in sorted(adata.obs['batch'].unique()):
        m = batch_vals == b
        ax.scatter(xy[m, 0], xy[m, 1], c=BATCH_COLORS[b],
                   s=0.3, alpha=0.4, linewidths=0, rasterized=True)
    leg = [Patch(color=BATCH_COLORS[b], label=b)
           for b in sorted(adata.obs['batch'].unique())]
    ax.legend(handles=leg, fontsize=6.5, frameon=False, ncol=1,
              loc='lower right', handlelength=0.8, labelspacing=0.3)
    ax.set_title(f'{method} — Dataset', fontsize=10, pad=3)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.set_xlabel('UMAP 1', fontsize=8); ax.set_ylabel('UMAP 2', fontsize=8)
    ax.tick_params(labelsize=6.5)
    ax.text(-0.07, 1.03, 'A' if row == 0 else 'D',
            transform=ax.transAxes, fontsize=13, fontweight='bold', va='bottom')

    # Disease
    ax = axes[row, 1]
    for dis in ['UC', 'CRC']:
        m = dis_vals == dis
        ax.scatter(xy[m, 0], xy[m, 1], c=DISEASE_COLORS[dis],
                   s=0.3, alpha=0.45, linewidths=0, rasterized=True)
    leg = [Patch(color=DISEASE_COLORS[d], label=d) for d in ['CRC', 'UC']]
    ax.legend(handles=leg, fontsize=8, frameon=False, loc='lower right',
              handlelength=0.8)
    ax.set_title(f'{method} — Disease', fontsize=10, pad=3)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.set_xlabel('UMAP 1', fontsize=8); ax.set_ylabel('UMAP 2', fontsize=8)
    ax.tick_params(labelsize=6.5)
    ax.text(-0.07, 1.03, 'B' if row == 0 else 'E',
            transform=ax.transAxes, fontsize=13, fontweight='bold', va='bottom')

    # Cell type
    ax = axes[row, 2]
    for ct in CT_ORDER:
        m = ct_vals == ct
        if m.sum() == 0:
            continue
        ax.scatter(xy[m, 0], xy[m, 1], c=CT_COLORS.get(ct, '#aaa'),
                   s=0.3, alpha=0.45, linewidths=0, rasterized=True)
    leg = [Patch(color=CT_COLORS[ct], label=ct) for ct in CT_ORDER
           if ct in adata.obs['cellass_type'].values]
    ax.legend(handles=leg, fontsize=6, frameon=False, ncol=1,
              loc='lower right', handlelength=0.8, labelspacing=0.25)
    ax.set_title(f'{method} — Cell Type', fontsize=10, pad=3)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.set_xlabel('UMAP 1', fontsize=8); ax.set_ylabel('UMAP 2', fontsize=8)
    ax.tick_params(labelsize=6.5)
    ax.text(-0.07, 1.03, 'C' if row == 0 else 'F',
            transform=ax.transAxes, fontsize=13, fontweight='bold', va='bottom')

fig.suptitle('scVI vs scANVI UMAP Comparison', fontsize=12, y=0.97)

for ext in ['png', 'pdf']:
    p = f'{FIG_DIR}/fig_scanvi_umap.{ext}'
    fig.savefig(p, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f"  Saved {FIG_DIR}/fig_scanvi_umap.png/.pdf")

print(f"\n[{datetime.now():%H:%M:%S}] scANVI training complete.")
print(f"  X_scvi  → X_umap_scvi   (original scVI UMAP)")
print(f"  X_scanvi → X_umap_scanvi (scANVI UMAP)")
print(f"  Both embeddings in: {out_h5ad}")
