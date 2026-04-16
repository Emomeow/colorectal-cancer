#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplementary QC figures: nCount, nFeature, MT% distributions for all datasets.
Phase 1 (5 main datasets) + Phase 2 (GSE144735) + Phase 3 (GSE178341).
Shows raw distribution (gray) overlaid with post-filter distribution (colored),
with threshold lines where applicable.

Output: figures/phase1/suppfig_qc_all_datasets.png / .pdf
        Project_Delivery/results/Supplementary/Supp1_QC_Metrics/figures/
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
os.makedirs('Project_Delivery/results/Supplementary/Supp1_QC_Metrics/figures', exist_ok=True)

# ── Dataset configs ────────────────────────────────────────────────────────────
# phase: 1 = Phase 1 main integration; 2 = Phase 2 spatial; 3 = Phase 3 MIL
DATASETS = [
    # ── Phase 1 ──
    dict(name='GSE132257', raw='processed_data/GSE132257.h5ad',
         pre='processed_data/GSE132257_preprocessed.h5ad',
         color='#e06c75', disease='CRC', phase=1,
         thresholds=dict(mt=10)),
    dict(name='GSE132465', raw='processed_data/GSE132465.h5ad',
         pre='processed_data/GSE132465_preprocessed.h5ad',
         color='#c678dd', disease='CRC', phase=1,
         thresholds=dict(mt=10)),
    dict(name='GSE125527', raw='processed_data/GSE125527.h5ad',
         pre='processed_data/GSE125527_preprocessed.h5ad',
         color='#61afef', disease='UC', phase=1,
         thresholds=dict(mt=10)),
    dict(name='GSE150115', raw='processed_data/GSE150115.h5ad',
         pre='processed_data/GSE150115_preprocessed.h5ad',
         color='#56b6c2', disease='UC', phase=1,
         thresholds=dict(mt=10, min_genes=500, min_counts=1000)),
    dict(name='cellxgene', raw='processed_data/cellxgene_raw.h5ad',
         pre='processed_data/cellxgene_preprocessed.h5ad',
         color='#e5c07b', disease='CRC', phase=1,
         thresholds=dict(mt=10, min_genes=200, min_counts=500)),
    # ── Phase 2 ──
    # GSE144735: paper-provided pre-filtered; use same file for raw/pre
    dict(name='GSE144735', raw='processed_data/GSE144735.h5ad',
         pre='processed_data/GSE144735.h5ad',
         color='#2CA02C', disease='CRC (spatial)', phase=2,
         note='paper-provided',
         thresholds=dict(mt=10)),
    # ── Phase 3 ──
    # GSE178341: MMRp h5ad (QC computed, pre-MT-filter) → immune_mil (MT<25%)
    dict(name='GSE178341', raw='processed_data/GSE178341_MMRp.h5ad',
         pre='processed_data/GSE178341_immune_mil.h5ad',
         color='#FF7F0E', disease='CRC (MMR)', phase=3,
         thresholds=dict(mt=25, min_genes=200)),
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

# Phase separator rows (drawn as horizontal lines between phase groups)
PHASE_LABELS = {1: 'Phase 1 — Main Integration (5 datasets)',
                2: 'Phase 2 — Spatial Analysis',
                3: 'Phase 3 — MIL Training'}

fig, axes = plt.subplots(N_DS, N_M, figsize=(14, 3.0 * N_DS),
                         gridspec_kw=dict(hspace=0.60, wspace=0.38))

prev_phase = None
for row_i, (cfg, raw_obs, pre_obs) in enumerate(records):
    color      = cfg['color']
    name       = cfg['name']
    thresholds = cfg['thresholds']
    phase      = cfg['phase']
    note       = cfg.get('note', '')
    same_file  = cfg['raw'] == cfg['pre']   # GSE144735: paper-provided

    for col_i, (metric, label, scale, thresh_key) in enumerate(METRICS):
        ax = axes[row_i, col_i]

        raw_vals = raw_obs[metric].values
        pre_vals = pre_obs[metric].values

        # Bin edges on raw data
        if scale == 'log':
            lo   = np.log10(max(raw_vals.min(), 1))
            hi   = np.log10(raw_vals.max()) + 0.05
            bins = np.logspace(lo, hi, 50)
        else:
            bins = np.linspace(0, raw_vals.max() * 1.05, 50)

        # Raw histogram (gray) — skip for paper-provided single-file datasets
        if not same_file:
            ax.hist(raw_vals, bins=bins, color='#cccccc', alpha=0.7,
                    density=True, linewidth=0)
        # Post-filter histogram (colored)
        ax.hist(pre_vals, bins=bins, color=color, alpha=0.85,
                density=True, linewidth=0)

        if scale == 'log':
            ax.set_xscale('log')

        # Threshold lines
        if thresh_key == 'mt' and 'mt' in thresholds:
            tv = thresholds['mt']
            ax.axvline(tv, color='#e74c3c', lw=1.4, ls='--', alpha=0.9)
            ax.text(tv * 1.03, ax.get_ylim()[1] * 0.85 if ax.get_ylim()[1] else 1,
                    f'{tv}%', color='#e74c3c', fontsize=7, va='top')
        elif thresh_key == 'min_counts' and 'min_counts' in thresholds:
            ax.axvline(thresholds['min_counts'], color='#e74c3c', lw=1.4, ls='--', alpha=0.9)
        elif thresh_key == 'min_genes' and 'min_genes' in thresholds:
            ax.axvline(thresholds['min_genes'], color='#e74c3c', lw=1.4, ls='--', alpha=0.9)

        # Median annotation
        med = np.median(pre_vals)
        med_str = (f'med={med:,.0f}' if scale == 'log'
                   else (f'med={med:.1f}%' if 'mt' in metric else f'med={med:.0f}'))
        ax.text(0.98, 0.97, med_str, transform=ax.transAxes,
                ha='right', va='top', fontsize=7.5, color='#2c3e50')

        # Column header (row 0 only)
        if row_i == 0:
            ax.set_title(label, fontsize=10, fontweight='bold', pad=6)

        # Y-axis label with phase tag
        if col_i == 0:
            n_raw = len(raw_vals)
            n_pre = len(pre_vals)
            if same_file:
                count_str = f'n={n_pre:,} (pre-filtered)'
            else:
                count_str = f'n={n_raw:,}→{n_pre:,}'
            phase_str = f'[Ph{phase}]'
            ax.set_ylabel(
                f"{phase_str} {name}\n({cfg['disease']})\n{count_str}",
                fontsize=7.5, labelpad=4
            )
            # Draw phase separator line above first row of each new phase
            if phase != prev_phase and row_i > 0:
                ax.axhline(y=ax.get_ylim()[1] if ax.get_ylim()[1] else 1,
                           color='#999999', lw=0, ls='-')  # placeholder; use spines
                # Add a visible top border to separate phases
                for spine in ax.spines.values():
                    pass
                ax.spines['top'].set_visible(True)
                ax.spines['top'].set_linewidth(1.5)
                ax.spines['top'].set_color('#888888')
                ax.spines['top'].set_linestyle('--')
                # Phase label as text above the row
                ax.text(-0.18, 1.18, PHASE_LABELS[phase],
                        transform=ax.transAxes, fontsize=8,
                        fontweight='bold', color='#555555',
                        va='bottom', ha='left')
        else:
            ax.set_ylabel('')
        ax.set_xlabel('')
        ax.tick_params(labelsize=7)
        ax.spines['right'].set_visible(False)
        if not (phase != prev_phase and row_i > 0 and col_i == 0):
            ax.spines['top'].set_visible(False)
        ax.relim(); ax.autoscale_view()

    prev_phase = phase

# ── Legend ─────────────────────────────────────────────────────────────────────
legend_handles = [
    mpatches.Patch(color='#cccccc', alpha=0.7, label='Before filtering'),
    mpatches.Patch(color='#888888', alpha=0.85, label='After filtering'),
    Line2D([0], [0], color='#e74c3c', lw=1.4, ls='--', label='Filter threshold'),
    mpatches.Patch(color='white', label='[Ph1]=Phase 1  [Ph2]=Phase 2  [Ph3]=Phase 3'),
]
fig.legend(handles=legend_handles, loc='upper center', ncol=4,
           fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, 1.01),
           bbox_transform=fig.transFigure)

fig.suptitle('Supplementary Figure S1 — QC Metrics Across All Datasets (Phases 1–3)',
             fontsize=12, fontweight='bold', y=1.04)

# ── Save ──────────────────────────────────────────────────────────────────────
OUTS = [
    'figures/phase1/suppfig_qc_all_datasets',
    'Project_Delivery/results/Supplementary/Supp1_QC_Metrics/figures/suppfig_qc_all_datasets',
]
for base in OUTS:
    for ext in ['png', 'pdf']:
        fig.savefig(f'{base}.{ext}', dpi=200, bbox_inches='tight',
                    format=ext if ext == 'pdf' else None)
        print(f"Saved: {base}.{ext}")

plt.close(fig)
print("Done.")
