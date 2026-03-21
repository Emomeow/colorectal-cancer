#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot Panel E & F using pre-computed scib-metrics results.
Re-computes bootstrap AMI (fast) for Panel F.
"""

import scanpy as sc
import numpy as np
import pandas as pd
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

# ── Hardcoded scib-metrics results (from completed benchmark run) ──────────────
# Individual metric scores per embedding
METRIC_DATA = {
    #                        scVI      Harmony   PCA(uninteg.)
    'Isolated labels':      [0.5806,   0.5677,   0.5844],
    'KMeans NMI':           [0.6388,   0.6114,   0.6518],
    'KMeans ARI':           [0.4776,   0.4233,   0.4420],
    'Silhouette label':     [0.5098,   0.5567,   0.5538],
    'cLISI':                [0.9947,   0.9943,   0.9953],
    'iLISI':                [0.0136,   0.0869,   0.0000],
    'kBET':                 [0.1002,   0.3273,   0.0627],
    'Graph connectivity':   [0.8653,   0.8495,   0.8388],
    'PCR comparison':       [0.8424,   0.4320,   0.0000],
}
# Which metrics belong to which category
BIO_METRICS   = ['cLISI', 'Silhouette label', 'KMeans NMI', 'KMeans ARI', 'Isolated labels']
BATCH_METRICS = ['iLISI', 'kBET', 'Graph connectivity', 'PCR comparison']
METRIC_ORDER  = BATCH_METRICS + BIO_METRICS   # batch on top, bio below

METHOD_ORDER  = ['scVI', 'Harmony', 'PCA\n(uninteg.)']
DISPLAY_ORDER = ['X_scvi', 'X_harmony', 'X_pca_uninteg']
NAME_MAP = {
    'X_scvi':        'scVI',
    'X_harmony':     'Harmony',
    'X_pca_uninteg': 'PCA\n(uninteg.)',
}
AGGREGATE = {
    'scVI':         {'Bio': 0.640, 'Batch': 0.455},
    'Harmony':      {'Bio': 0.631, 'Batch': 0.424},
    'PCA\n(uninteg.)': {'Bio': 0.645, 'Batch': 0.225},
}

# ── Bootstrap cluster stability (fast, ~2 min) ────────────────────────────────
N_BOOT      = 10
RESOLUTIONS = [0.3, 0.5, 0.8, 1.0]

print(f"[{datetime.now():%H:%M:%S}] Loading data for bootstrap...")
adata = sc.read_h5ad('scvi_results_hvg3k/merged_scvi_hvg3k.h5ad')
print(f"  {adata.n_obs:,} cells")

print(f"[{datetime.now():%H:%M:%S}] Bootstrap stability...")
if 'neighbors' not in adata.uns:
    sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=15)

boot_results       = {}
n_clusters_per_res = {}
for res in RESOLUTIONS:
    base_key = f'leiden_r{res:.1f}'
    sc.tl.leiden(adata, resolution=res, key_added=base_key, random_state=0)
    base_labs = adata.obs[base_key].values.astype(str)
    n_cl = adata.obs[base_key].nunique()
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

# ── Panel E: heatmap of all scib-metrics ─────────────────────────────────────
ax = axes[0]

# Build matrix: rows=metrics (METRIC_ORDER), cols=methods
mat = np.array([[METRIC_DATA[m][i] for i in range(3)] for m in METRIC_ORDER])

import matplotlib.colors as mcolors
cmap = plt.cm.RdYlGn
im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=1, aspect='auto')

# Cell value annotations
for i in range(len(METRIC_ORDER)):
    for j in range(3):
        v = mat[i, j]
        txt_color = 'white' if (v < 0.25 or v > 0.80) else 'black'
        ax.text(j, i, f'{v:.3f}', ha='center', va='center',
                fontsize=7, color=txt_color, fontweight='normal')

# Separator line between batch and bio metrics
sep = len(BATCH_METRICS) - 0.5
ax.axhline(sep, color='white', lw=2)

# Category labels on left margin
ax.text(-0.52, (len(BATCH_METRICS) - 1) / 2, 'Batch\ncorrection',
        transform=ax.get_yaxis_transform(), ha='center', va='center',
        fontsize=7, color='#3C5488', fontweight='bold', rotation=90)
ax.text(-0.52, len(BATCH_METRICS) + (len(BIO_METRICS) - 1) / 2, 'Bio\nconservation',
        transform=ax.get_yaxis_transform(), ha='center', va='center',
        fontsize=7, color='#E64B35', fontweight='bold', rotation=90)

# Aggregate score row label
for j, name in enumerate(METHOD_ORDER):
    bio_v   = AGGREGATE[name]['Bio']
    batch_v = AGGREGATE[name]['Batch']
    ax.text(j, len(METRIC_ORDER) + 0.05,
            f'Bio={bio_v:.3f}\nBatch={batch_v:.3f}',
            ha='center', va='bottom', fontsize=6, color='#444')

ax.set_xticks(range(3))
ax.set_xticklabels(METHOD_ORDER, fontsize=8.5)
ax.set_yticks(range(len(METRIC_ORDER)))
ax.set_yticklabels(METRIC_ORDER, fontsize=7.5)
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')
ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False, labelsize=8)

# Colorbar
cbar = fig.colorbar(im, ax=ax, shrink=0.55, aspect=15, pad=0.02)
cbar.set_label('Score (0–1)', fontsize=7)
cbar.ax.tick_params(labelsize=6.5)

ax.text(-0.12, 1.03, 'E', transform=ax.transAxes,
        fontsize=14, fontweight='bold', va='bottom')
ax.set_title('Integration Quality (scib-metrics)', fontsize=10, pad=22)

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
ax.text(-0.10, 1.03, 'F', transform=ax.transAxes,
        fontsize=14, fontweight='bold', va='bottom')
ax.set_title('Cluster Stability (Bootstrap AMI)', fontsize=10, pad=4)

for ext in ['png', 'pdf']:
    fig.savefig(f'figures/phase1/fig1_panelEF.{ext}',
                dpi=200, bbox_inches='tight')
plt.close(fig)
print(f"\n[{datetime.now():%H:%M:%S}] Saved figures/phase1/fig1_panelEF.png/.pdf")
