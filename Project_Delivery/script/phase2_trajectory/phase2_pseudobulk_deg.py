#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 2 — Pseudobulk DEG (CRC vs UC) with pyDESeq2 + GSEA
Per-subtype: aggregate raw counts → DESeq2 → volcano → GSEA (Hallmark)
"""

import scanpy as sc
import scvi
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.sparse import issparse
import gseapy as gp
from statsmodels.stats.multitest import multipletests
import os, sys, io, warnings
from datetime import datetime

warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

plt.rcParams['font.family']  = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

FIG_DIR = 'figures/phase2/pseudobulk'
RES_DIR = 'phase2_results/pseudobulk'
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

BATCH_UC  = {'GSE125527', 'GSE150115'}

# Subtypes to test (biologically prioritised)
TARGET_SUBTYPES = [
    'CD8_Tex', 'TAM', 'Treg', 'CD8_Teff', 'CD8_Tpex',
    'Classical_Mono', 'Th17', 'cDC2', 'CD4_Tconv',
    'NK cells', 'B cells',
]
MIN_CELLS_PER_SAMPLE = 10
MIN_SAMPLES_PER_COND = 5   # relaxed to 3 for TAM
MIN_SAMPLES_TAM      = 3


def savefig(fig, name):
    fig.savefig(f'{FIG_DIR}/{name}.png', dpi=200, bbox_inches='tight', facecolor='white')
    fig.savefig(f'{FIG_DIR}/{name}.pdf',           bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'   Saved: {FIG_DIR}/{name}.png/.pdf')


# ══════════════════════════════════════════════════════════════════
#  scVI batch-aware differential expression
# ══════════════════════════════════════════════════════════════════

ALL_BATCHES = ['GSE132465', 'GSE132257', 'GSE125527', 'GSE150115', 'cellxgene']
MAX_CELLS_PER_GROUP = 1500   # cap for memory efficiency


def run_scvi_de(model, adata, subtype, min_cells=100):
    """
    Batch-corrected DE using scVI get_normalized_expression(transform_batch).

    For each cell, the expression is decoded through ALL batches and averaged,
    effectively marginalising out batch effects. Then Wilcoxon rank-sum test
    gives gene-level p-values and log2FC is computed from the corrected means.
    """
    from scipy.stats import ranksums

    obs = adata.obs
    mask_crc = (obs['cell_subtype_fine'] == subtype) & (obs['condition'] == 'CRC')
    mask_uc  = (obs['cell_subtype_fine'] == subtype) & (obs['condition'] == 'UC')

    n_crc, n_uc = mask_crc.sum(), mask_uc.sum()
    if n_crc < min_cells or n_uc < min_cells:
        return None

    # Cap cell numbers for memory/speed
    idx_crc = np.where(mask_crc)[0]
    idx_uc  = np.where(mask_uc)[0]
    rng = np.random.default_rng(42)
    if len(idx_crc) > MAX_CELLS_PER_GROUP:
        idx_crc = rng.choice(idx_crc, MAX_CELLS_PER_GROUP, replace=False)
    if len(idx_uc) > MAX_CELLS_PER_GROUP:
        idx_uc  = rng.choice(idx_uc,  MAX_CELLS_PER_GROUP, replace=False)

    print(f'   scVI DE: {len(idx_crc):,} CRC cells vs {len(idx_uc):,} UC cells  (batch-corrected)')

    # get_normalized_expression with transform_batch marginalises over all batches
    adata_crc = adata[idx_crc]
    adata_uc  = adata[idx_uc]

    expr_crc = model.get_normalized_expression(
        adata           = adata_crc,
        transform_batch = ALL_BATCHES,
        n_samples       = 10,
        return_numpy    = True,
    )   # shape: (n_cells, n_genes) — averaged over samples
    expr_uc = model.get_normalized_expression(
        adata           = adata_uc,
        transform_batch = ALL_BATCHES,
        n_samples       = 10,
        return_numpy    = True,
    )

    genes    = adata.var_names.tolist()
    mean_crc = expr_crc.mean(axis=0)
    mean_uc  = expr_uc.mean(axis=0)
    log2fc   = np.log2(mean_crc + 1e-6) - np.log2(mean_uc + 1e-6)

    # Wilcoxon rank-sum per gene (vectorised over genes is slow; use only genes with any signal)
    pvals = np.ones(len(genes))
    for j in range(expr_crc.shape[1]):
        c, u = expr_crc[:, j], expr_uc[:, j]
        if c.max() > 0 or u.max() > 0:
            _, pvals[j] = ranksums(c, u)

    _, fdr, _, _ = multipletests(pvals, method='fdr_bh')

    df = pd.DataFrame({
        'log2FC':   log2fc,
        'pvalue':   pvals,
        'fdr':      fdr,
        'mean_crc': mean_crc,
        'mean_uc':  mean_uc,
    }, index=genes)
    df.index.name = 'gene'
    df = df.sort_values('fdr')
    return df


# ══════════════════════════════════════════════════════════════════
#  Volcano plot
# ══════════════════════════════════════════════════════════════════

def plot_volcano(res, title, name, fc_cut=1.0, fdr_cut=0.05, top_n=20):
    fig, ax = plt.subplots(figsize=(7, 6))
    fig.suptitle(title, fontsize=11, fontweight='bold')

    log2fc = res['log2FC'].values
    neglog = -np.log10(res['fdr'].clip(1e-300).values)
    sig_up   = (res['fdr'] < fdr_cut) & (res['log2FC'] >  fc_cut)
    sig_down = (res['fdr'] < fdr_cut) & (res['log2FC'] < -fc_cut)
    ns       = ~(sig_up | sig_down)

    ax.scatter(log2fc[ns],       neglog[ns],       c='#CCCCCC', s=4,  alpha=0.4, linewidths=0)
    ax.scatter(log2fc[sig_up],   neglog[sig_up],   c='#E64B35', s=8,  alpha=0.7, linewidths=0, label=f'CRC↑ ({sig_up.sum()})')
    ax.scatter(log2fc[sig_down], neglog[sig_down], c='#4DBBD5', s=8,  alpha=0.7, linewidths=0, label=f'UC↑ ({sig_down.sum()})')

    ax.axvline( fc_cut, color='grey', lw=0.8, ls='--')
    ax.axvline(-fc_cut, color='grey', lw=0.8, ls='--')
    ax.axhline(-np.log10(fdr_cut), color='grey', lw=0.8, ls='--')
    ax.set_xlabel('log₂FC (CRC / UC)', fontsize=10)
    ax.set_ylabel('-log₁₀(FDR)', fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=8)

    # Label top genes
    top = pd.concat([
        res[sig_up].head(top_n // 2),
        res[sig_down].head(top_n // 2)
    ])
    for gene, row in top.iterrows():
        ax.text(row['log2FC'], -np.log10(row['fdr'].clip(1e-300)) + 0.3,
                gene, fontsize=6, ha='center', va='bottom',
                color='#B22222' if row['log2FC'] > 0 else '#1A6B8A')

    plt.tight_layout()
    savefig(fig, name)

    return sig_up.sum(), sig_down.sum()


# ══════════════════════════════════════════════════════════════════
#  GSEA (Hallmark gene sets via gseapy)
# ══════════════════════════════════════════════════════════════════

def run_gsea(res, subtype_tag, top_n_plot=15):
    """Pre-ranked GSEA on Hallmark gene sets."""
    # Rank genes by sign(log2FC) * -log10(pvalue)
    rank_metric = res['log2FC'] * (-np.log10(res['pvalue'].clip(1e-300)))
    rank_metric = rank_metric.sort_values(ascending=False)
    rank_metric.name = 'score'

    try:
        gsea_res = gp.prerank(
            rnk          = rank_metric,
            gene_sets    = 'MSigDB_Hallmark_2020',
            outdir       = None,
            min_size     = 10,
            max_size     = 500,
            permutation_num = 500,
            seed         = 42,
            verbose      = False,
        )
        df = gsea_res.res2d.copy()
        df = df.rename(columns={'NES': 'NES', 'FDR q-val': 'fdr', 'Name': 'pathway'})
        if 'pathway' not in df.columns and df.index.name == 'Term':
            df = df.reset_index().rename(columns={'Term': 'pathway'})
        df['fdr'] = pd.to_numeric(df['fdr'], errors='coerce')
        df['NES'] = pd.to_numeric(df['NES'],  errors='coerce')
        df = df.dropna(subset=['NES', 'fdr']).sort_values('fdr')
        df.to_csv(f'{RES_DIR}/gsea_{subtype_tag}.csv', index=False)

        # Plot top pathways
        plot_gsea_bar(df, subtype_tag, top_n=top_n_plot)
        return df
    except Exception as e:
        print(f'   GSEA failed: {e}')
        return None


def plot_gsea_bar(df, subtype_tag, top_n=15):
    """Horizontal bar plot of top GSEA hits."""
    sig = df[df['fdr'] < 0.25].copy()
    if sig.empty:
        sig = df.head(top_n)

    # Top positive + negative NES
    top_pos = sig[sig['NES'] > 0].head(top_n // 2)
    top_neg = sig[sig['NES'] < 0].head(top_n // 2)
    plot_df  = pd.concat([top_pos, top_neg]).sort_values('NES')

    if plot_df.empty:
        return

    fig, ax = plt.subplots(figsize=(8, max(4, len(plot_df) * 0.42)))
    fig.suptitle(f'GSEA Hallmark — {subtype_tag}', fontsize=11, fontweight='bold')

    colors = ['#E64B35' if v > 0 else '#4DBBD5' for v in plot_df['NES']]
    bars = ax.barh(range(len(plot_df)), plot_df['NES'].values,
                   color=colors, alpha=0.85, edgecolor='none')

    path_labels = plot_df['pathway'].str.replace('HALLMARK_', '', regex=False).str.replace('_', ' ', regex=False)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(path_labels.tolist(), fontsize=8)
    ax.axvline(0, color='black', lw=0.8)
    ax.set_xlabel('NES (Normalized Enrichment Score)', fontsize=9)
    ax.spines[['top', 'right']].set_visible(False)

    # FDR annotations
    for i, (_, row) in enumerate(plot_df.iterrows()):
        fdr_str = f'q={row["fdr"]:.3f}' if row['fdr'] >= 0.001 else f'q={row["fdr"]:.1e}'
        ax.text(row['NES'] + (0.05 if row['NES'] > 0 else -0.05),
                i, fdr_str, fontsize=6, va='center',
                ha='left' if row['NES'] > 0 else 'right')

    plt.tight_layout()
    savefig(fig, f'gsea_{subtype_tag}')


# ══════════════════════════════════════════════════════════════════
#  Summary heatmap: top DEGs across subtypes
# ══════════════════════════════════════════════════════════════════

def plot_deg_summary_heatmap(all_results, name, top_n=10):
    """Heatmap of log2FC for top DEGs across all subtypes."""
    # Collect top up + down genes per subtype
    gene_set = set()
    for st, res in all_results.items():
        sig = res[(res['fdr'] < 0.05) & (res['log2FC'].abs() > 1)]
        top_up   = sig[sig['log2FC'] > 0].head(top_n).index.tolist()
        top_down = sig[sig['log2FC'] < 0].head(top_n).index.tolist()
        gene_set.update(top_up + top_down)

    genes = sorted(gene_set)
    if not genes:
        return

    subtypes = list(all_results.keys())
    mat = pd.DataFrame(index=genes, columns=subtypes, dtype=float)
    for st, res in all_results.items():
        for g in genes:
            if g in res.index:
                mat.loc[g, st] = res.loc[g, 'log2FC']

    mat = mat.fillna(0)
    # Sort genes by mean log2FC
    mat = mat.loc[mat.mean(axis=1).sort_values(ascending=False).index]

    vmax = np.percentile(mat.abs().values, 95)
    fig, ax = plt.subplots(figsize=(max(6, len(subtypes) * 1.1),
                                     max(6, len(genes) * 0.28)))
    fig.suptitle('Top DEGs across subtypes (log₂FC CRC/UC)', fontsize=11, fontweight='bold')

    im = ax.imshow(mat.values, aspect='auto', cmap='RdBu_r',
                   vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(subtypes)))
    ax.set_xticklabels(subtypes, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=7)
    plt.colorbar(im, ax=ax, shrink=0.6, label='log₂FC')
    plt.tight_layout()
    savefig(fig, name)


# ══════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════
print('=' * 65)
print(f'Phase 2 — scVI DE  |  {datetime.now().strftime("%H:%M:%S")}')
print('=' * 65)

# ── Load data ────────────────────────────────────────────────
print('\nLoading data ...')
adata = sc.read_h5ad('processed_data/merged_annotated_fine.h5ad')
print(f'  {adata.n_obs:,} cells')

# Add condition column (CRC / UC) based on batch
adata.obs['condition'] = adata.obs['batch'].apply(
    lambda b: 'UC' if b in BATCH_UC else 'CRC'
)
print(f"  CRC: {(adata.obs['condition']=='CRC').sum():,}  UC: {(adata.obs['condition']=='UC').sum():,}")

# ── Load scVI model ──────────────────────────────────────────
print('\nLoading scVI model ...')
model = scvi.model.SCVI.load('scvi_model_optimized/', adata=adata)
print('  Model loaded.')

# ── Per-subtype DE ───────────────────────────────────────────
all_results  = {}
all_n_sig_up = {}
all_n_sig_dn = {}

MIN_CELLS     = 100   # minimum cells per condition for scVI DE
MIN_CELLS_TAM =  30   # relaxed for TAM (CRC-enriched, few UC cells)

for subtype in TARGET_SUBTYPES:
    print(f'\n{"─"*55}')
    print(f'[{subtype}]  {datetime.now().strftime("%H:%M:%S")}')

    n_sub = (adata.obs['cell_subtype_fine'] == subtype).sum()
    if n_sub < 50:
        print(f'  Skipped: only {n_sub} cells')
        continue

    try:
        mc = MIN_CELLS_TAM if subtype == 'TAM' else MIN_CELLS
        res = run_scvi_de(model, adata, subtype, min_cells=mc)
        if res is None:
            print(f'  Skipped: insufficient cells per condition')
            continue

        tag = subtype.replace(' ', '_')
        res.to_csv(f'{RES_DIR}/deg_{tag}.csv')
        all_results[subtype] = res

        n_up = ((res['fdr'] < 0.05) & (res['log2FC'] > 1)).sum()
        n_dn = ((res['fdr'] < 0.05) & (res['log2FC'] < -1)).sum()
        all_n_sig_up[subtype] = n_up
        all_n_sig_dn[subtype] = n_dn
        print(f'   DEGs: {n_up} CRC↑, {n_dn} UC↑  (FDR<0.05, |log2FC|>1)')

        # Top genes
        top_crc = res[(res['fdr'] < 0.05) & (res['log2FC'] > 1)].head(10)
        top_uc  = res[(res['fdr'] < 0.05) & (res['log2FC'] < -1)].head(10)
        if not top_crc.empty:
            print(f'   CRC↑ top: {", ".join(top_crc.index[:8].tolist())}')
        if not top_uc.empty:
            print(f'   UC↑  top: {", ".join(top_uc.index[:8].tolist())}')

        # Volcano
        plot_volcano(res, f'DEG: {subtype}  (CRC vs UC)', f'volcano_{tag}')

        # GSEA
        print(f'   Running GSEA (Hallmark) ...')
        run_gsea(res, tag)

    except Exception as e:
        print(f'   Failed: {e}')
        continue

# ── Summary ──────────────────────────────────────────────────
print('\n' + '=' * 65)
print('SUMMARY')
print('=' * 65)
summary = pd.DataFrame({
    'CRC_up': all_n_sig_up,
    'UC_up':  all_n_sig_dn,
}).fillna(0).astype(int)
print(summary.to_string())
summary.to_csv(f'{RES_DIR}/deg_summary.csv')

if len(all_results) >= 2:
    print('\nPlotting cross-subtype heatmap ...')
    plot_deg_summary_heatmap(all_results, 'deg_heatmap_all_subtypes', top_n=8)

print(f'\nDone  |  {datetime.now().strftime("%H:%M:%S")}')
print('=' * 65)
