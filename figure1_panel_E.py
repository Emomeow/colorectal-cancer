#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 1 Panel E — Integration Quality: scVI vs Harmony vs Unintegrated (PCA)
  iLISI : batch mixing (higher = better)
  NMI   : clustering concordance with known cell types (higher = better)
  ARI   : clustering concordance with known cell types (higher = better)
All metrics computed on immune cells only (T/B/NK/Myeloids/Mast cells).
"""
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys, io
from datetime import datetime

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

os.makedirs('figures/phase1', exist_ok=True)

print("=" * 60)
print(f"Figure 1E — Integration Benchmark  |  {datetime.now().strftime('%H:%M:%S')}")
print("=" * 60)

# ─────────────────────────────────────────────────────────────────
# 1. Load data
# ─────────────────────────────────────────────────────────────────
print("\n[1/5] Loading data ...")
adata = sc.read_h5ad('processed_data/merged_annotated.h5ad')
print(f"  {adata.n_obs:,} cells  |  X_scvi: {adata.obsm['X_scvi'].shape}")

batch_labels    = adata.obs['batch'].astype(str).values
celltype_labels = adata.obs['cell_type_grouped'].astype(str).values

# ─────────────────────────────────────────────────────────────────
# 2. Build PCA embedding (unintegrated baseline)
# ─────────────────────────────────────────────────────────────────
print("\n[2/5] Computing PCA (unintegrated baseline) ...")
adata_pca = adata.copy()
sc.pp.normalize_total(adata_pca, target_sum=1e4)
sc.pp.log1p(adata_pca)
sc.pp.highly_variable_genes(adata_pca, n_top_genes=3000, batch_key='batch')
sc.pp.scale(adata_pca, max_value=10)
sc.tl.pca(adata_pca, n_comps=30, use_highly_variable=True)
X_pca = adata_pca.obsm['X_pca']
print(f"  PCA shape: {X_pca.shape}")

# ─────────────────────────────────────────────────────────────────
# 3. Build Harmony embedding (call harmonypy directly; scanpy wrapper
#    has a shape bug with harmonypy >= 0.2.0)
# ─────────────────────────────────────────────────────────────────
print("\n[3/5] Running Harmony ...")
import harmonypy as hm
ho = hm.run_harmony(
    X_pca.astype(np.float64),
    pd.DataFrame({'batch': batch_labels}),
    'batch',
    max_iter_harmony=20,
    random_state=0,
    verbose=True,
)
X_harmony = np.array(ho.Z_corr)   # shape (n_cells, 30)
print(f"  Harmony shape: {X_harmony.shape}")

# ─────────────────────────────────────────────────────────────────
# 4. Subset to immune cells
#
#  Why immune cells only?
#  cellxgene contains ONLY malignant + epithelial cells (no immune).
#  GSE125527 contains ALMOST NO malignant/epithelial cells.
#  iLISI on the full dataset is dominated by these non-overlapping
#  populations whose neighbourhood is 100 % single-batch.
#  Immune cells are present in all 5 batches and are the correct
#  population to evaluate cross-batch mixing.
# ─────────────────────────────────────────────────────────────────
print("\n[4/5] Computing metrics (immune cells only) ...")
from scib_metrics.nearest_neighbors import pynndescent
import scib_metrics

IMMUNE_TYPES = {'T cells', 'NK cells', 'B cells', 'Myeloids', 'Mast cells'}
immune_mask  = adata.obs['cell_type_grouped'].isin(IMMUNE_TYPES).values
print(f"  Immune cells: {immune_mask.sum():,} / {len(immune_mask):,} total")
print("  Batch breakdown (immune cells):")
for b, n in adata[immune_mask].obs['batch'].value_counts().items():
    print(f"    {b}: {n:,}")

batch_labels_imm    = batch_labels[immune_mask]
celltype_labels_imm = celltype_labels[immune_mask]

N_NEIGHBORS = 90

results = {}
for name, X_full in [('Unintegrated\n(PCA)', X_pca),
                      ('Harmony',             X_harmony),
                      ('scVI',                adata.obsm['X_scvi'])]:
    X_emb = X_full[immune_mask].astype(np.float32)
    print(f"\n  [{name.strip()}]", flush=True)

    nn_res = pynndescent(X_emb, n_neighbors=N_NEIGHBORS, random_state=0, n_jobs=-1)

    # iLISI — batch mixing (higher = better)
    ilisi = scib_metrics.ilisi_knn(nn_res, batch_labels_imm, scale=True)
    print(f"    iLISI      = {ilisi:.4f}")

    # NMI / ARI — Leiden clustering concordance with known cell types
    nmi_ari = scib_metrics.nmi_ari_cluster_labels_leiden(nn_res, celltype_labels_imm)
    # nmi_ari is a dict with keys 'nmi' and 'ari'
    nmi = float(nmi_ari['nmi'])
    ari = float(nmi_ari['ari'])
    print(f"    NMI        = {nmi:.4f}")
    print(f"    ARI        = {ari:.4f}")

    results[name] = {'iLISI': ilisi, 'NMI': nmi, 'ARI': ari}

# ─────────────────────────────────────────────────────────────────
# 5. Plot: 3-panel bar chart (iLISI | NMI | ARI)
# ─────────────────────────────────────────────────────────────────
print("\n[5/5] Plotting ...")

METHOD_ORDER  = ['Unintegrated\n(PCA)', 'Harmony', 'scVI']
METHOD_LABELS = ['Unintegrated\n(PCA)', 'Harmony', 'scVI']
METHOD_COLORS = ['#AAAAAA', '#4DBBD5', '#E64B35']

n_immune = immune_mask.sum()
fig, axes = plt.subplots(1, 3, figsize=(12, 4.5))
fig.suptitle(
    f'Figure 1E — Integration Quality (immune cells, n={n_immune:,})',
    fontsize=12, fontweight='bold', y=1.03
)

PANELS = [
    ('iLISI',  'iLISI — Batch mixing\n(higher is better)',        'Scale 0–1'),
    ('NMI',    'NMI — Cell type concordance\n(higher is better)', 'NMI'),
    ('ARI',    'ARI — Cell type concordance\n(higher is better)', 'ARI'),
]

x = np.arange(len(METHOD_ORDER))

for ax, (metric, title, ylabel) in zip(axes, PANELS):
    vals = [results[m][metric] for m in METHOD_ORDER]
    bars = ax.bar(x, vals, width=0.55, color=METHOD_COLORS,
                  edgecolor='white', linewidth=0.6)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.008,
                f'{v:.3f}', ha='center', va='bottom',
                fontsize=8.5, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(METHOD_LABELS, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_ylim(0, max(vals) * 1.30)
    ax.set_title(title, fontsize=10, fontweight='bold', pad=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='y', labelsize=8)

plt.tight_layout()

out_png = 'figures/phase1/figure1_panel_E_lisi.png'
out_pdf = 'figures/phase1/figure1_panel_E_lisi.pdf'
plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(out_pdf, bbox_inches='tight', facecolor='white')
plt.close()
print(f"\nSaved: {out_png}")
print(f"Saved: {out_pdf}")
print("Done.")
