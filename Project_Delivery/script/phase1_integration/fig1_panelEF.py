#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 1 – Panel E & F
  E: Integration quality via scib-metrics Benchmarker
     (scVI / Harmony / PCA-unintegrated, batch + bio conservation)
  F: Cluster stability: Bootstrap AMI at multiple Leiden resolutions
"""

import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_mutual_info_score
import warnings, os, sys, io
from datetime import datetime

warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
plt.rcParams['font.family'] = 'Arial'
os.makedirs('figures/phase1', exist_ok=True)

# ── Config ────────────────────────────────────────────────────────────────────
N_BOOT      = 10
RESOLUTIONS = [0.3, 0.5, 0.8, 1.0]

# ── Load ──────────────────────────────────────────────────────────────────────
print(f"[{datetime.now():%H:%M:%S}] Loading data...")
adata = sc.read_h5ad('scvi_results_hvg3k/merged_scvi_hvg3k.h5ad')
print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

# ── PCA baseline ──────────────────────────────────────────────────────────────
print(f"[{datetime.now():%H:%M:%S}] PCA on normalized data...")
adata_norm = ad.AnnData(
    X=adata.layers['counts'].copy(),
    obs=adata.obs.copy(),
    var=adata.var.copy(),
)
sc.pp.normalize_total(adata_norm, target_sum=1e4)
sc.pp.log1p(adata_norm)
sc.pp.scale(adata_norm, max_value=10)
sc.pp.pca(adata_norm, n_comps=30, svd_solver='arpack')
adata.obsm['X_pca_uninteg'] = adata_norm.obsm['X_pca'].copy()
print(f"  PCA: {adata.obsm['X_pca_uninteg'].shape}")

# ── Harmony ───────────────────────────────────────────────────────────────────
print(f"[{datetime.now():%H:%M:%S}] Running Harmony...")
try:
    import harmonypy as hm
    ho = hm.run_harmony(adata_norm.obsm['X_pca'], adata_norm.obs, 'batch',
                        max_iter_harmony=20, verbose=False)
    Z = np.array(ho.Z_corr)
    adata.obsm['X_harmony'] = Z.T if Z.shape[0] != adata.n_obs else Z
    has_harmony = True
    print(f"  Harmony: {adata.obsm['X_harmony'].shape}")
except Exception as e:
    print(f"  Harmony failed: {e}")
    has_harmony = False

# ── scib-metrics Benchmarker ──────────────────────────────────────────────────
print(f"[{datetime.now():%H:%M:%S}] Running scib-metrics Benchmarker...")
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

embedding_keys = ['X_pca_uninteg', 'X_scvi']
if has_harmony:
    embedding_keys.insert(1, 'X_harmony')  # order: PCA, Harmony, scVI

bio_metrics = BioConservation(
    isolated_labels=True,
    nmi_ari_cluster_labels_leiden=False,
    nmi_ari_cluster_labels_kmeans=True,
    silhouette_label=True,
    clisi_knn=True,
)
batch_metrics = BatchCorrection(
    bras=False,
    ilisi_knn=True,
    kbet_per_label=True,
    graph_connectivity=True,
    pcr_comparison=True,
)

bm = Benchmarker(
    adata,
    batch_key='batch',
    label_key='cellass_type',
    embedding_obsm_keys=embedding_keys,
    bio_conservation_metrics=bio_metrics,
    batch_correction_metrics=batch_metrics,
    pre_integrated_embedding_obsm_key='X_pca_uninteg',
    n_jobs=-1,
)
bm.benchmark()
results_df = bm.get_results(min_max_scale=False)
print("\nscib-metrics results:")
print(results_df.to_string())

# ── Extract aggregate bio conservation & batch correction ─────────────────────
# scib-metrics uses "Batch correction" and "Bio conservation" as aggregate cols
# Try several possible column names across versions
bio_cols   = [c for c in results_df.columns if 'bio' in c.lower()]
batch_cols = [c for c in results_df.columns if 'batch' in c.lower()]
print(f"\nBio columns: {bio_cols}")
print(f"Batch columns: {batch_cols}")

# Use the last (aggregate) bio and batch column
bio_col   = bio_cols[-1]   if bio_cols   else None
batch_col = batch_cols[-1] if batch_cols else None

# Map embedding key → display name
NAME_MAP = {
    'X_scvi':        'scVI',
    'X_harmony':     'Harmony',
    'X_pca_uninteg': 'PCA\n(uninteg.)',
}
# Display order: scVI first (best expected)
DISPLAY_ORDER = ['X_scvi', 'X_harmony', 'X_pca_uninteg']
DISPLAY_ORDER = [k for k in DISPLAY_ORDER if k in results_df.index]

print(f"\nBio col: {bio_col}  |  Batch col: {batch_col}")
for k in DISPLAY_ORDER:
    bio_v   = results_df.loc[k, bio_col]   if bio_col   else float('nan')
    batch_v = results_df.loc[k, batch_col] if batch_col else float('nan')
    print(f"  {NAME_MAP[k]:20s}  bio={bio_v:.3f}  batch={batch_v:.3f}")

# ── Bootstrap cluster stability ───────────────────────────────────────────────
print(f"\n[{datetime.now():%H:%M:%S}] Bootstrap stability "
      f"(res={RESOLUTIONS}, n_boot={N_BOOT})...")
if 'neighbors' not in adata.uns:
    sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=15)

boot_results       = {}
n_clusters_per_res = {}
for res in RESOLUTIONS:
    base_key = f'leiden_r{res:.1f}'
    sc.tl.leiden(adata, resolution=res, key_added=base_key, random_state=0)
    base_labs = adata.obs[base_key].values.astype(str)
    n_cl      = adata.obs[base_key].nunique()
    n_clusters_per_res[res] = n_cl
    amis = []
    for seed in range(N_BOOT):
        sc.tl.leiden(adata, resolution=res, key_added='_boot',
                     random_state=(seed + 1) * 13)
        ami = adjusted_mutual_info_score(base_labs,
                                         adata.obs['_boot'].values.astype(str))
        amis.append(ami)
    adata.obs.drop(columns=['_boot'], inplace=True, errors='ignore')
    boot_results[res] = np.array(amis)
    print(f"  res={res:.1f}  n_clusters={n_cl:3d}  "
          f"AMI={boot_results[res].mean():.4f} ± {boot_results[res].std():.4f}")

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
fig.subplots_adjust(left=0.07, right=0.97, top=0.90, bottom=0.13, wspace=0.40)

# ── Panel E ───────────────────────────────────────────────────────────────────
ax = axes[0]

if bio_col and batch_col:
    methods   = DISPLAY_ORDER
    labels    = [NAME_MAP[k] for k in methods]
    bio_vals  = [results_df.loc[k, bio_col]   for k in methods]
    batch_vals= [results_df.loc[k, batch_col] for k in methods]

    x = np.arange(len(methods))
    W = 0.34
    b1 = ax.bar(x - W/2, batch_vals, W, color='#4DBBD5', alpha=0.88,
                label='Batch correction (↑ better)')
    b2 = ax.bar(x + W/2, bio_vals,   W, color='#E64B35', alpha=0.88,
                label='Bio conservation (↑ better)')

    for bar, v in zip(b1, batch_vals):
        ax.text(bar.get_x() + bar.get_width()/2, v + 0.008,
                f'{v:.3f}', ha='center', va='bottom', fontsize=6.5)
    for bar, v in zip(b2, bio_vals):
        ax.text(bar.get_x() + bar.get_width()/2, v + 0.008,
                f'{v:.3f}', ha='center', va='bottom', fontsize=6.5)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ymax = max(max(bio_vals), max(batch_vals)) + 0.18
    ax.set_ylim(0, min(ymax, 1.22))
else:
    ax.text(0.5, 0.5, 'scib-metrics results\nnot available',
            transform=ax.transAxes, ha='center', va='center', fontsize=9)

ax.set_ylabel('Aggregate score (scib-metrics)', fontsize=9)
ax.axhline(1.0, color='#ccc', lw=0.8, ls='--', alpha=0.6)
ax.tick_params(labelsize=7.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(fontsize=7.5, frameon=False, loc='upper left',
          handlelength=0.8, labelspacing=0.4)
ax.text(-0.08, 1.03, 'F', transform=ax.transAxes,
        fontsize=14, fontweight='bold', va='bottom')
ax.set_title('Integration Quality (scib-metrics)', fontsize=10, pad=4)
ax.text(0.98, 0.02, f'n={adata.n_obs:,} cells\n10 CellAssign types',
        transform=ax.transAxes, ha='right', va='bottom',
        fontsize=6.5, color='#666')

# ── Panel F ───────────────────────────────────────────────────────────────────
ax = axes[1]
x_pos = np.arange(len(RESOLUTIONS))
rng2  = np.random.default_rng(123)

bp = ax.boxplot(
    [boot_results[r] for r in RESOLUTIONS],
    positions=x_pos,
    widths=0.45,
    patch_artist=True,
    medianprops=dict(color='#E64B35', linewidth=2.2),
    whiskerprops=dict(color='#555', linewidth=1),
    capprops=dict(color='#555', linewidth=1),
    flierprops=dict(marker='', markersize=0),
    boxprops=dict(facecolor='#91D1C2', alpha=0.75, linewidth=0.8),
)
for xi, res in zip(x_pos, RESOLUTIONS):
    jitter = rng2.normal(0, 0.07, N_BOOT)
    ax.scatter(xi + jitter, boot_results[res], s=38, c='#3C5488',
               alpha=0.88, zorder=5, linewidths=0.5, edgecolors='white')

xlabels = [f'res={r:.1f}\n({n_clusters_per_res[r]} clusters)'
           for r in RESOLUTIONS]
ax.set_xticks(x_pos)
ax.set_xticklabels(xlabels, fontsize=7.5)
ax.set_xlim(-0.6, len(RESOLUTIONS) - 0.4)
ami_min = min(v.min() for v in boot_results.values())
ax.set_ylim(max(0, ami_min - 0.03), 1.04)
ax.set_ylabel('Adjusted Mutual Information (AMI)', fontsize=9)
ax.tick_params(labelsize=7.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.axhline(1.0, color='#ccc', lw=0.8, ls='--', alpha=0.6)
ax.text(0.97, 0.05,
        f'n={N_BOOT} bootstrap seeds\nscVI Leiden clustering',
        transform=ax.transAxes, ha='right', va='bottom',
        fontsize=7, color='#444')
ax.text(-0.10, 1.03, '', transform=ax.transAxes,
        fontsize=14, fontweight='bold', va='bottom')
ax.set_title('Cluster Stability (Bootstrap AMI)', fontsize=10, pad=4)

for ext in ['png', 'pdf']:
    fig.savefig(f'figures/phase1/fig1_panelEF.{ext}',
                dpi=200, bbox_inches='tight')
plt.close(fig)
print(f"\n[{datetime.now():%H:%M:%S}] Saved figures/phase1/fig1_panelEF.png/.pdf")
