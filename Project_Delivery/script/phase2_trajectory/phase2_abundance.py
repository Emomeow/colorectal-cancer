#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 2 — Differential Cell Type Abundance (CRC vs UC)
Proportion-based analysis: per-sample subtype fractions → Mann-Whitney + BH FDR
Uses fine annotations from phase2_fine_annotation.py
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import os, sys, io, warnings, re
from datetime import datetime

warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

plt.rcParams['font.family']  = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

FIG_DIR = 'figures/phase2/abundance'
RES_DIR = 'phase2_results/abundance'
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

BATCH_UC = {'GSE125527', 'GSE150115'}


def savefig(fig, name):
    fig.savefig(f'{FIG_DIR}/{name}.png', dpi=200, bbox_inches='tight', facecolor='white')
    fig.savefig(f'{FIG_DIR}/{name}.pdf',           bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'   Saved: {FIG_DIR}/{name}.png/.pdf')


# ══════════════════════════════════════════════════════════════════
#  1. Build sample_id column
# ══════════════════════════════════════════════════════════════════

def build_sample_ids(adata):
    """
    Extract patient/sample-level IDs from existing metadata and barcodes.
    Use numpy arrays throughout to avoid pandas index-alignment issues.
    """
    n = adata.n_obs
    sample_id_arr = np.full(n, '', dtype=object)
    condition_arr = np.full(n, 'exclude', dtype=object)

    batch_arr = adata.obs['batch'].values  # numpy array (no index)
    barcodes  = np.array(adata.obs_names.tolist())

    # ── GSE125527 ─────────────────────────────────────────────
    mask = batch_arr == 'GSE125527'
    bc_series = pd.Series(barcodes[mask])
    prefix = bc_series.str.extract(r'^([A-Z]\d+)_', expand=False).fillna('unk').values
    sample_id_arr[mask] = np.array(['GSE125527_' + p for p in prefix])
    condition_arr[mask] = 'UC'

    # ── GSE150115 ─────────────────────────────────────────────
    mask = batch_arr == 'GSE150115'
    bc_series = pd.Series(barcodes[mask])
    gsm = bc_series.str.extract(r'^(GSM\d+)_', expand=False).fillna('unk').values
    sample_id_arr[mask] = gsm
    condition_arr[mask] = 'UC'

    # ── GSE132465 ─────────────────────────────────────────────
    mask = batch_arr == 'GSE132465'
    patient = adata.obs['Patient'].values[mask].astype(str)
    sample_id_arr[mask] = np.array(['GSE132465_' + p for p in patient])
    condition_arr[mask] = 'CRC'

    # ── GSE132257 ─────────────────────────────────────────────
    mask = batch_arr == 'GSE132257'
    patient = adata.obs['Patient'].values[mask].astype(str)
    sample_id_arr[mask] = np.array(['GSE132257_' + p for p in patient])
    condition_arr[mask] = 'CRC'

    adata.obs['sample_id'] = sample_id_arr
    adata.obs['condition'] = condition_arr
    return adata


# ══════════════════════════════════════════════════════════════════
#  2. Compute per-sample subtype proportions
# ══════════════════════════════════════════════════════════════════

def compute_proportions(adata, cell_type_col='cell_subtype_fine',
                         exclude_subtypes=None, min_cells=50):
    """
    For each sample, compute proportion of each subtype among immune cells.
    Returns:
      prop_df  — samples × subtypes (proportions)
      meta_df  — sample metadata (condition, batch, n_cells)
    """
    if exclude_subtypes is None:
        exclude_subtypes = {'Unknown', 'Unknown_T', 'Unknown_My',
                            'Epithelial', 'Malignant', 'Stromal'}

    keep = (adata.obs['condition'].isin(['CRC', 'UC'])) & \
           (~adata.obs[cell_type_col].isin(exclude_subtypes))
    sub  = adata[keep]

    # Count cells per sample × subtype
    ct = pd.crosstab(sub.obs['sample_id'], sub.obs[cell_type_col])
    ct = ct[ct.sum(axis=1) >= min_cells]   # drop tiny samples

    prop = ct.div(ct.sum(axis=1), axis=0)  # row-normalize → proportions

    # Sample metadata
    meta = sub.obs[['sample_id', 'condition', 'batch']].drop_duplicates('sample_id').set_index('sample_id')
    meta = meta.loc[prop.index]

    return prop, meta, ct


# ══════════════════════════════════════════════════════════════════
#  3. Mann-Whitney test per subtype
# ══════════════════════════════════════════════════════════════════

def test_proportions(prop, meta):
    """Mann-Whitney U test: CRC vs UC per subtype. Returns sorted result df."""
    crc_idx = meta[meta['condition'] == 'CRC'].index
    uc_idx  = meta[meta['condition'] == 'UC'].index

    rows = []
    for st in prop.columns:
        crc_vals = prop.loc[crc_idx, st].values
        uc_vals  = prop.loc[uc_idx,  st].values
        stat, pval = mannwhitneyu(crc_vals, uc_vals, alternative='two-sided')
        log2fc = np.log2((crc_vals.mean() + 1e-6) / (uc_vals.mean() + 1e-6))
        rows.append({
            'subtype': st,
            'CRC_mean': crc_vals.mean(),
            'UC_mean':  uc_vals.mean(),
            'log2FC':   log2fc,
            'pval':     pval,
            'n_CRC':    len(crc_vals),
            'n_UC':     len(uc_vals),
        })

    df = pd.DataFrame(rows)
    _, fdr, _, _ = multipletests(df['pval'], method='fdr_bh')
    df['fdr'] = fdr
    df = df.sort_values('fdr')
    return df


# ══════════════════════════════════════════════════════════════════
#  4. Plotting
# ══════════════════════════════════════════════════════════════════

SUBTYPE_COLORS = {
    # T cells
    'CD8_Teff':  '#E64B35', 'CD8_Tex':   '#B22222', 'CD8_Tpex':  '#FF8C69',
    'CD4_Tconv': '#4DBBD5', 'Treg':      '#00A087', 'Th17':      '#3C5488',
    'Tfh':       '#F39B7F', 'gdT':       '#8491B4', 'Cycling_T': '#91D1C2',
    # Myeloid
    'Classical_Mono': '#E64B35', 'NonClassical_Mono': '#F39B7F',
    'TAM':        '#7B0000', 'Inflam_Mac': '#FF6600',
    'cDC1':       '#4DBBD5', 'cDC2':      '#00A087',
    'pDC':        '#3C5488', 'Proliferating_My': '#91D1C2',
    # Other
    'B cells': '#AEC7E8', 'NK cells': '#FFBB78', 'Mast cells': '#C5B0D5',
}


def plot_dotplot_abundance(result_df, title, name):
    """
    Bubble plot: x = log2FC(CRC/UC), y = subtype (sorted by log2FC)
    Bubble size = -log10(FDR), color = direction
    """
    df = result_df.copy().sort_values('log2FC')
    df['-log10FDR'] = -np.log10(df['fdr'].clip(1e-10))
    sig = df['fdr'] < 0.05

    fig, ax = plt.subplots(figsize=(7, max(5, len(df) * 0.4)))
    fig.suptitle(title, fontsize=12, fontweight='bold')

    if df.empty:
        ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center')
        plt.tight_layout()
        savefig(fig, name)
        return

    for i, row in df.reset_index(drop=True).iterrows():
        col = '#E64B35' if row['log2FC'] > 0 else '#4DBBD5'
        alpha = 0.9 if row['fdr'] < 0.05 else 0.4
        size  = row['-log10FDR'] * 80 + 30
        ax.scatter(row['log2FC'], i, s=size, c=col, alpha=alpha,
                   linewidths=0.5, edgecolors='k' if row['fdr'] < 0.05 else 'grey',
                   zorder=3)

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['subtype'].tolist(), fontsize=9)
    ax.axvline(0, color='grey', lw=0.8, ls='--', zorder=1)
    ax.set_xlabel('log₂FC (CRC / UC)', fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(axis='x', color='#EEEEEE', lw=0.5, zorder=0)

    # Significance markers
    for i, row in df.reset_index(drop=True).iterrows():
        if row['fdr'] < 0.05:
            ax.text(row['log2FC'] + 0.05, i, '★', fontsize=8,
                    color='#E64B35' if row['log2FC'] > 0 else '#4DBBD5',
                    va='center')

    # Size legend
    for fdr_v, lbl in [(0.05, 'FDR=0.05'), (0.01, 'FDR=0.01'), (0.001, 'FDR=0.001')]:
        s = -np.log10(fdr_v) * 80 + 30
        ax.scatter([], [], s=s, c='grey', alpha=0.7, linewidths=0, label=lbl)
    ax.legend(title='Bubble size', fontsize=7, title_fontsize=7,
              loc='lower right')

    plt.tight_layout()
    savefig(fig, name)


def plot_stripplot(prop, meta, result_df, cell_group, name, top_n=12):
    """
    Strip plot (jitter) + box: proportion per sample, split CRC vs UC.
    Show top_n most significant subtypes.
    """
    sig_subtypes = result_df.sort_values('fdr').head(top_n)['subtype'].tolist()
    # keep only those present in prop
    sig_subtypes = [s for s in sig_subtypes if s in prop.columns]
    if not sig_subtypes:
        return

    n = len(sig_subtypes)
    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.2, nrows * 3.2))
    fig.suptitle(f'Differential abundance — {cell_group} (top {n} by FDR)',
                 fontsize=11, fontweight='bold')
    axes = np.array(axes).flatten()

    rng = np.random.default_rng(42)
    for ax, st in zip(axes, sig_subtypes):
        row = result_df[result_df['subtype'] == st].iloc[0]
        for cond, col, x in [('UC', '#4DBBD5', 0), ('CRC', '#E64B35', 1)]:
            idx  = meta[meta['condition'] == cond].index
            vals = prop.loc[idx, st].values * 100
            jit  = rng.uniform(-0.15, 0.15, len(vals))
            ax.scatter([x + j for j in jit], vals, c=col, s=25, alpha=0.7,
                       linewidths=0, zorder=3)
            ax.plot([x - 0.2, x + 0.2], [np.median(vals)] * 2,
                    color='black', lw=2, zorder=4)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['UC', 'CRC'], fontsize=9)
        ax.set_ylabel('Proportion (%)', fontsize=8)
        fdr_str = f'FDR={row["fdr"]:.3f}' if row['fdr'] >= 0.001 else f'FDR={row["fdr"]:.1e}'
        fc_str  = f'log₂FC={row["log2FC"]:+.2f}'
        ax.set_title(f'{st}\n{fc_str}  {fdr_str}', fontsize=8, fontweight='bold')
        ax.spines[['top', 'right']].set_visible(False)

    for ax in axes[n:]:
        ax.axis('off')

    plt.tight_layout()
    savefig(fig, name)


def plot_composition_bar(prop, meta, title, name):
    """Stacked bar: mean proportion per subtype, CRC vs UC side by side."""
    crc_mean = prop.loc[meta[meta['condition'] == 'CRC'].index].mean()
    uc_mean  = prop.loc[meta[meta['condition'] == 'UC'].index].mean()

    # Sort subtypes by CRC_mean - UC_mean (most CRC-enriched first)
    order = (crc_mean - uc_mean).sort_values(ascending=False).index.tolist()
    crc_vals = crc_mean[order]
    uc_vals  = uc_mean[order]

    fig, ax = plt.subplots(figsize=(max(8, len(order) * 0.55), 5))
    fig.suptitle(title, fontsize=11, fontweight='bold')

    x = np.arange(len(order))
    w = 0.35
    for i, st in enumerate(order):
        col = SUBTYPE_COLORS.get(st, '#999999')
        ax.bar(x[i] - w/2, crc_vals[st] * 100, w, color=col, alpha=0.9, label=st if i == 0 else '')
        ax.bar(x[i] + w/2, uc_vals[st]  * 100, w, color=col, alpha=0.45)

    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Mean proportion (%)', fontsize=9)

    from matplotlib.patches import Patch
    legend_els = [Patch(facecolor='grey', alpha=0.9, label='CRC'),
                  Patch(facecolor='grey', alpha=0.4, label='UC')]
    ax.legend(handles=legend_els, fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    savefig(fig, name)


# ══════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════
print('=' * 65)
print(f'Phase 2 — Differential Abundance  |  {datetime.now().strftime("%H:%M:%S")}')
print('=' * 65)

print('\nLoading fine-annotated data ...')
adata = sc.read_h5ad('processed_data/merged_annotated_fine.h5ad')
print(f'  {adata.n_obs:,} cells')

# ── Build sample IDs ──────────────────────────────────────────
print('\n[1] Building sample IDs ...')
adata = build_sample_ids(adata)

sample_summary = adata.obs[adata.obs['condition'].isin(['CRC','UC'])].groupby(
    ['condition', 'batch', 'sample_id']).size().reset_index(name='n_cells')
print(f'\n  Samples per condition:')
print(sample_summary.groupby('condition').agg(n_samples=('sample_id','nunique'),
                                               n_cells=('n_cells','sum')).to_string())
print(f'\n  CRC samples: {sample_summary[sample_summary.condition=="CRC"]["sample_id"].nunique()}')
print(f'  UC samples:  {sample_summary[sample_summary.condition=="UC"]["sample_id"].nunique()}')

# ── All immune cells (T + B + NK + Myeloid + Mast) ───────────
print('\n[2] All-immune-cell abundance ...')
immune_types = {'T cells', 'B cells', 'NK cells', 'Myeloids', 'Mast cells',
                'CD8_Teff', 'CD8_Tex', 'CD8_Tpex', 'CD4_Tconv', 'Treg',
                'Th17', 'Tfh', 'gdT', 'Cycling_T',
                'Classical_Mono', 'NonClassical_Mono', 'TAM', 'Inflam_Mac',
                'cDC1', 'cDC2', 'pDC', 'Proliferating_My'}

prop_all, meta_all, ct_all = compute_proportions(adata, min_cells=50)
print(f'  {len(prop_all)} samples × {len(prop_all.columns)} subtypes')
print(f'  CRC: {(meta_all.condition=="CRC").sum()} samples  |  UC: {(meta_all.condition=="UC").sum()} samples')

result_all = test_proportions(prop_all, meta_all)
result_all.to_csv(f'{RES_DIR}/abundance_all_subtypes.csv', index=False)

print('\n  Significant subtypes (FDR < 0.05):')
sig = result_all[result_all['fdr'] < 0.05]
for _, row in sig.iterrows():
    direction = 'CRC↑' if row['log2FC'] > 0 else 'UC↑'
    print(f'    {row["subtype"]:<22} log2FC={row["log2FC"]:+.2f}  FDR={row["fdr"]:.3e}  [{direction}]')

plot_dotplot_abundance(result_all, 'Differential abundance — all immune subtypes',
                       'abundance_dotplot_all')

# ── T cells only ──────────────────────────────────────────────
print('\n[3] T cell abundance ...')
t_subtypes = ['CD8_Teff', 'CD8_Tex', 'CD8_Tpex', 'CD4_Tconv',
              'Treg', 'Th17', 'Tfh', 'gdT', 'Cycling_T']

adata_t = adata[adata.obs['cell_type_grouped'] == 'T cells'].copy()
adata_t.obs['condition'] = adata.obs.loc[adata_t.obs_names, 'condition']
adata_t.obs['sample_id'] = adata.obs.loc[adata_t.obs_names, 'sample_id']

prop_t, meta_t, ct_t = compute_proportions(adata_t, cell_type_col='cell_subtype_fine',
                                             min_cells=30)
result_t = test_proportions(prop_t, meta_t)
result_t.to_csv(f'{RES_DIR}/abundance_tcell.csv', index=False)

print('\n  T cell results (sorted by FDR):')
print(result_t[['subtype','log2FC','pval','fdr','CRC_mean','UC_mean']].to_string(index=False))

plot_dotplot_abundance(result_t, 'Differential abundance — T cells', 'abundance_dotplot_tcell')
plot_stripplot(prop_t, meta_t, result_t, 'T cells', 'abundance_strip_tcell', top_n=9)
plot_composition_bar(prop_t, meta_t, 'T cell subtype composition CRC vs UC', 'abundance_bar_tcell')

# ── Myeloid cells only ────────────────────────────────────────
print('\n[4] Myeloid cell abundance ...')
my_subtypes = ['Classical_Mono', 'NonClassical_Mono', 'TAM', 'Inflam_Mac',
               'cDC1', 'cDC2', 'pDC', 'Proliferating_My']

adata_my = adata[adata.obs['cell_type_grouped'] == 'Myeloids'].copy()
adata_my.obs['condition'] = adata.obs.loc[adata_my.obs_names, 'condition']
adata_my.obs['sample_id'] = adata.obs.loc[adata_my.obs_names, 'sample_id']

prop_my, meta_my, ct_my = compute_proportions(adata_my, cell_type_col='cell_subtype_fine',
                                               min_cells=10)
result_my = test_proportions(prop_my, meta_my)
result_my.to_csv(f'{RES_DIR}/abundance_myeloid.csv', index=False)

print('\n  Myeloid results (sorted by FDR):')
print(result_my[['subtype','log2FC','pval','fdr','CRC_mean','UC_mean']].to_string(index=False))

plot_dotplot_abundance(result_my, 'Differential abundance — Myeloid', 'abundance_dotplot_myeloid')
plot_stripplot(prop_my, meta_my, result_my, 'Myeloid', 'abundance_strip_myeloid', top_n=8)
plot_composition_bar(prop_my, meta_my, 'Myeloid subtype composition CRC vs UC', 'abundance_bar_myeloid')

# ── Summary ──────────────────────────────────────────────────
print('\n' + '=' * 65)
print('SUMMARY — Significantly differential subtypes (FDR<0.05):')
print('=' * 65)
all_sig = result_all[result_all['fdr'] < 0.05].copy()
all_sig['direction'] = all_sig['log2FC'].apply(lambda x: 'CRC enriched' if x > 0 else 'UC enriched')
print(all_sig[['subtype','log2FC','fdr','direction']].sort_values('log2FC', ascending=False).to_string(index=False))

print(f'\nDone  |  {datetime.now().strftime("%H:%M:%S")}')
print('=' * 65)
