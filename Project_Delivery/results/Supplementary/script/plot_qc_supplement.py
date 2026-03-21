#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplementary QC figures: nCount, nFeature, MT% distributions for all 5 datasets.
Shows raw distribution (gray) overlaid with post-filter distribution (colored),
with threshold lines where applicable.

Output: figures/phase1/suppfig_qc_all_datasets.png / .pdf
"""

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import os

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
plt.rcParams['font.family'] = 'Arial'

os.makedirs('figures/phase1', exist_ok=True)

# ── Dataset configs ────────────────────────────────────────────────────────────
DATASETS = [
    dict(name='GSE132257', raw='processed_data/GSE132257.h5ad',
         pre='processed_data/GSE132257_preprocessed.h5ad',
         color='#e06c75', disease='CRC',
         thresholds=dict(mt=20)),
    dict(name='GSE132465', raw='processed_data/GSE132465.h5ad',
         pre='processed_data/GSE132465_preprocessed.h5ad',
         color='#c678dd', disease='CRC',
         thresholds=dict(mt=20)),
    dict(name='GSE125527', raw='processed_data/GSE125527.h5ad',
         pre='processed_data/GSE125527_preprocessed.h5ad',
         color='#61afef', disease='UC',
         thresholds=dict(mt=20)),
    dict(name='GSE150115', raw='processed_data/GSE150115.h5ad',
         pre='processed_data/GSE150115_preprocessed.h5ad',
         color='#56b6c2', disease='UC',
         thresholds=dict(mt=20, min_genes=500, min_counts=1000)),
    dict(name='cellxgene', raw='processed_data/cellxgene_raw.h5ad',
         pre='processed_data/cellxgene_preprocessed.h5ad',
         color='#e5c07b', disease='CRC',
         thresholds=dict(mt=20, min_genes=200, min_counts=500)),
]

METRICS = [
    ('total_counts',      'nCount (UMI)',     'log',  'min_counts'),
    ('n_genes_by_counts', 'nFeature (genes)', 'log',  'min_genes'),
    ('pct_counts_mt',     'MT%',              'linear', 'mt'),
]

# ── Load data ─────────────────────────────────────────────────────────────────
print("Loading datasets...")

def get_obs(path, calc_qc=False):
    adata = ad.read_h5ad(path)
    if 'n_genes_by_counts' not in adata.obs.columns:
        import scanpy as sc
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                                   log1p=False, inplace=True)
    return adata.obs[['total_counts', 'n_genes_by_counts', 'pct_counts_mt']].copy()

records = []
for cfg in DATASETS:
    print(f"  {cfg['name']}...")
    raw_obs  = get_obs(cfg['raw'])
    pre_obs  = get_obs(cfg['pre'])
    raw_obs['_stage'] = 'raw'
    pre_obs['_stage'] = 'clean'
    combined = pd.concat([raw_obs, pre_obs])
    combined['_dataset'] = cfg['name']
    records.append((cfg, raw_obs, pre_obs))

print("Plotting...")

# ── Figure layout ─────────────────────────────────────────────────────────────
N_DS = len(DATASETS)
N_M  = len(METRICS)

fig, axes = plt.subplots(N_DS, N_M, figsize=(14, 3.2 * N_DS),
                         gridspec_kw=dict(hspace=0.55, wspace=0.38))

for row_i, (cfg, raw_obs, pre_obs) in enumerate(records):
    color = cfg['color']
    name  = cfg['name']
    thresholds = cfg['thresholds']

    for col_i, (metric, label, scale, thresh_key) in enumerate(METRICS):
        ax = axes[row_i, col_i]

        raw_vals  = raw_obs[metric].values
        pre_vals  = pre_obs[metric].values

        # Compute shared bin edges on raw data
        if scale == 'log':
            lo = np.log10(max(raw_vals.min(), 1))
            hi = np.log10(raw_vals.max()) + 0.05
            bins = np.logspace(lo, hi, 50)
        else:
            lo = 0
            hi = raw_vals.max() * 1.05
            bins = np.linspace(lo, hi, 50)

        # Raw histogram (gray background)
        ax.hist(raw_vals,  bins=bins, color='#cccccc', alpha=0.7,
                label='Before', density=True, linewidth=0)
        # Post-filter histogram (colored foreground)
        ax.hist(pre_vals, bins=bins, color=color,     alpha=0.85,
                label='After',  density=True, linewidth=0)

        if scale == 'log':
            ax.set_xscale('log')

        # Threshold line
        thresh_val = None
        if thresh_key == 'mt' and 'mt' in thresholds:
            thresh_val = thresholds['mt']
            ax.axvline(thresh_val, color='#e74c3c', lw=1.4, ls='--', alpha=0.9)
            ax.text(thresh_val * 1.02, ax.get_ylim()[1] * 0.92 if ax.get_ylim()[1] > 0 else 1,
                    f'{thresh_val}%', color='#e74c3c', fontsize=7, va='top')
        elif thresh_key == 'min_counts' and 'min_counts' in thresholds:
            thresh_val = thresholds['min_counts']
            ax.axvline(thresh_val, color='#e74c3c', lw=1.4, ls='--', alpha=0.9)
        elif thresh_key == 'min_genes' and 'min_genes' in thresholds:
            thresh_val = thresholds['min_genes']
            ax.axvline(thresh_val, color='#e74c3c', lw=1.4, ls='--', alpha=0.9)

        # Statistics annotations
        n_raw  = len(raw_vals)
        n_pre  = len(pre_vals)
        med_pre = np.median(pre_vals)
        if scale == 'log':
            med_str = f'med={med_pre:,.0f}'
        else:
            med_str = f'med={med_pre:.1f}%' if 'mt' in metric.lower() else f'med={med_pre:.0f}'

        ax.text(0.98, 0.97, med_str,
                transform=ax.transAxes, ha='right', va='top',
                fontsize=7.5, color='#2c3e50')

        # Labels
        if row_i == 0:
            ax.set_title(label, fontsize=10, fontweight='bold', pad=6)
        if col_i == 0:
            disease_tag = cfg['disease']
            ax.set_ylabel(
                f"{name}\n({disease_tag})\nn={n_raw:,}→{n_pre:,}",
                fontsize=8, labelpad=4
            )
        else:
            ax.set_ylabel('')
        ax.set_xlabel('')

        ax.tick_params(labelsize=7)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Recalculate ylim after plotting so threshold text lands correctly
        ax.relim()
        ax.autoscale_view()

# ── Legend ─────────────────────────────────────────────────────────────────────
legend_handles = [
    mpatches.Patch(color='#cccccc', alpha=0.7, label='Before cleaning'),
    mpatches.Patch(color='#888888', alpha=0.85, label='After cleaning'),
    Line2D([0], [0], color='#e74c3c', lw=1.4, ls='--', label='Filter threshold'),
]
fig.legend(handles=legend_handles, loc='upper center', ncol=3,
           fontsize=9, frameon=False, bbox_to_anchor=(0.5, 1.0),
           bbox_transform=fig.transFigure)

fig.suptitle('Supplementary Figure S1 — QC Metrics: Before and After Cell Filtering',
             fontsize=12, fontweight='bold', y=1.02)

# ── Save ──────────────────────────────────────────────────────────────────────
for ext in ['png', 'pdf']:
    out = f'figures/phase1/suppfig_qc_all_datasets.{ext}'
    fig.savefig(out, dpi=200, bbox_inches='tight',
                format=ext if ext == 'pdf' else None)
    print(f"Saved: {out}")

plt.close(fig)
print("Done.")
