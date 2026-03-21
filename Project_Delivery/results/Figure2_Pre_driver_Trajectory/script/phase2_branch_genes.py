#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 2 — Branch Point Gene Analysis
Identify genes that drive T cell vs Myeloid fate divergence at the branch point.
Two complementary approaches:
  A. Palantir gene trends  — smooth expression curves along each lineage
  B. DEG at branch point   — Wilcoxon between fate-committed cells in the
                              high-entropy zone (the divergence window)
"""
import scanpy as sc
import numpy as np
import pandas as pd
import palantir
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys, io
from datetime import datetime

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

os.makedirs('figures/phase2', exist_ok=True)
os.makedirs('phase2_results', exist_ok=True)

print("=" * 65)
print(f"Phase 2 — Branch Point Gene Analysis  |  {datetime.now().strftime('%H:%M:%S')}")
print("=" * 65)

# ─────────────────────────────────────────────────────────────────
# 1. Load saved pseudotime adata + results
# ─────────────────────────────────────────────────────────────────
print("\n[1/5] Loading data ...")
adata = sc.read_h5ad('phase2_results/adata_pseudotime.h5ad')
pr_pseudotime   = pd.read_csv('phase2_results/pseudotime_values.csv',   index_col=0).squeeze()
pr_branch_probs = pd.read_csv('phase2_results/branch_probs.csv',         index_col=0)

# Map to standard Palantir obs keys (required by presults API)
adata.obs['palantir_pseudotime'] = pr_pseudotime.reindex(adata.obs_names).values
adata.obsm['palantir_fate_probabilities'] = pr_branch_probs.reindex(adata.obs_names)

# Terminal labels: T cell branch and Myeloid branch
T_TERM  = pr_branch_probs.columns[0]   # T cells terminal
MY_TERM = pr_branch_probs.columns[1]   # Myeloids terminal
print(f"  T cell terminal : {T_TERM}")
print(f"  Myeloid terminal: {MY_TERM}")

# ─────────────────────────────────────────────────────────────────
# 2. Normalize expression + select HVGs
# ─────────────────────────────────────────────────────────────────
print("\n[2/5] Normalizing expression and selecting HVGs ...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVGs computed on the full subset (batch-aware)
sc.pp.highly_variable_genes(adata, n_top_genes=2500, batch_key='batch')
adata_hvg = adata[:, adata.var['highly_variable']].copy()
print(f"  HVG subset: {adata_hvg.n_obs:,} cells × {adata_hvg.n_vars:,} genes")

# Copy Palantir keys to HVG subset
adata_hvg.obs['palantir_pseudotime'] = adata.obs['palantir_pseudotime']
adata_hvg.obsm['palantir_fate_probabilities'] = adata.obsm['palantir_fate_probabilities']
adata_hvg.obs['entropy']             = adata.obs['entropy']
adata_hvg.obsm['X_umap']             = adata.obsm['X_umap']

# ─────────────────────────────────────────────────────────────────
# 3A. Palantir gene trends along each lineage
# ─────────────────────────────────────────────────────────────────
print("\n[3/5] Computing gene trends (Palantir) ...")

# Select branch-committed cells (q=0.01 → bottom 1% fate-prob cells excluded)
palantir.presults.select_branch_cells(adata_hvg, q=0.01, eps=0.01)
print(f"  branch_masks stored in adata_hvg.uns")

# Compute gene trends (GAM spline fit per lineage)
# n_jobs=-1, limit to HVGs; this takes a few minutes
trends = palantir.presults.compute_gene_trends(
    adata_hvg,
    lineages=[T_TERM, MY_TERM],
    pseudo_time_key='palantir_pseudotime',
)
n_genes_with_trends = len(list(trends.values())[0]['trends'].columns) if trends else 0
print(f"  Gene trends computed: {n_genes_with_trends} genes × {len(trends)} lineages")

# ─────────────────────────────────────────────────────────────────
# 3B. Branch-point DEG: high-entropy zone, T-committed vs My-committed
# ─────────────────────────────────────────────────────────────────
print("\n[4/5] DEG: fate-committed cells (T-fate vs Myeloid-fate) ...")

# Use fate probability threshold to define committed cells:
# T-committed  : T_TERM prob > 0.7, any pseudotime
# My-committed : MY_TERM prob > 0.7, any pseudotime
# This gives much more cells than the high-entropy zone approach
pt_vals = adata.obs['palantir_pseudotime'].values
prob_T  = pr_branch_probs[T_TERM].reindex(adata.obs_names).values
prob_My = pr_branch_probs[MY_TERM].reindex(adata.obs_names).values

# Top 30% T-committed + all Myeloid-committed to keep groups balanced
T_THRESH  = np.percentile(prob_T,  70)   # top 30% T probability
MY_THRESH = np.percentile(prob_My, 50)   # top 50% My probability (all 107 + more)
t_committed  = prob_T  > T_THRESH
my_committed = prob_My > MY_THRESH
print(f"  T-committed  (prob>{T_THRESH:.2f}): {t_committed.sum():,} cells")
print(f"  My-committed (prob>{MY_THRESH:.2f}): {my_committed.sum():,} cells")

# Build DEG subset
keep = t_committed | my_committed
adata_bp = adata[keep].copy()
adata_bp.obs['fate_label'] = np.where(
    prob_T[keep] > prob_My[keep], 'T_cell_fate', 'Myeloid_fate'
)
n_T  = (adata_bp.obs['fate_label'] == 'T_cell_fate').sum()
n_My = (adata_bp.obs['fate_label'] == 'Myeloid_fate').sum()
print(f"  DEG subset: T_fate={n_T:,}  Myeloid_fate={n_My:,}")

# Wilcoxon does NOT need a KNN graph — remove stale uns keys to avoid conflicts
for stale_key in ['neighbors', 'diffmap_evals']:
    adata_bp.uns.pop(stale_key, None)

sc.tl.rank_genes_groups(
    adata_bp,
    groupby='fate_label',
    groups=['T_cell_fate', 'Myeloid_fate'],
    reference='rest',
    method='wilcoxon',
    key_added='branch_deg',
)
# Extract top DEGs
deg_T  = sc.get.rank_genes_groups_df(adata_bp, group='T_cell_fate',  key='branch_deg')
deg_My = sc.get.rank_genes_groups_df(adata_bp, group='Myeloid_fate', key='branch_deg')
deg_T  = deg_T[deg_T['logfoldchanges'] > 0].head(30)
deg_My = deg_My[deg_My['logfoldchanges'] > 0].head(30)

deg_T.to_csv('phase2_results/branch_DEG_Tcell_fate.csv', index=False)
deg_My.to_csv('phase2_results/branch_DEG_Myeloid_fate.csv', index=False)
print("  Top T-cell fate genes:", ', '.join(deg_T['names'].head(10).tolist()))
print("  Top Myeloid fate genes:", ', '.join(deg_My['names'].head(10).tolist()))

# ─────────────────────────────────────────────────────────────────
# 5. Plots
# ─────────────────────────────────────────────────────────────────
print("\n[5/5] Plotting ...")

# ── A: Gene trend heatmaps via Palantir ──────────────────────────
# Get genes that were actually modelled (may be subset of HVGs)
trend_genes = set()
for term_key, trend_dict in adata_hvg.uns.get('gene_trends', {}).items():
    if isinstance(trend_dict, dict) and 'trends' in trend_dict:
        trend_genes.update(trend_dict['trends'].columns.tolist())
print(f"  Genes with computed trends: {len(trend_genes)}")

try:
    palantir.presults.cluster_gene_trends(adata_hvg, branch_name=T_TERM,
                                          gene_trend_key='gene_trends')
    fig_hm = palantir.plot.plot_gene_trend_heatmaps(
        adata_hvg, gene_trend_key='gene_trends', scaling='z-score',
    )
    fig_hm.savefig('figures/phase2/gene_trend_heatmaps.png',
                   dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig_hm)
    print("  Saved: figures/phase2/gene_trend_heatmaps.png")
except Exception as e:
    print(f"  [WARN] Heatmap skipped: {e}")

# ── B: Top branch-point gene expression trends ───────────────────
# Use genes that overlap between DEG results and trend_genes
top_T_trend  = [g for g in deg_T['names'].head(12)  if g in trend_genes][:6]
top_My_trend = [g for g in deg_My['names'].head(12) if g in trend_genes][:6]
trend_plot_genes = top_T_trend + top_My_trend

if trend_plot_genes:
    try:
        fig_tr = palantir.plot.plot_gene_trends(
            adata_hvg, trend_plot_genes, gene_trend_key='gene_trends',
        )
        fig_tr.savefig('figures/phase2/branch_gene_trends.png',
                       dpi=200, bbox_inches='tight', facecolor='white')
        plt.close(fig_tr)
        print(f"  Saved: figures/phase2/branch_gene_trends.png  ({len(trend_plot_genes)} genes)")
    except Exception as e:
        print(f"  [WARN] Gene trends plot skipped: {e}")
else:
    print("  [WARN] No overlapping genes for trend plot")

# ── C: Branch-point DEG volcano plots ────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle('Branch-point DEG: T-cell fate vs Myeloid fate\n'
             '(cells in high-entropy branch zone)', fontsize=11, fontweight='bold')

for ax, df, title, color in [
    (ax1, deg_T,  'T-cell fate genes',  '#E64B35'),
    (ax2, deg_My, 'Myeloid fate genes', '#00A087'),
]:
    df_all = sc.get.rank_genes_groups_df(
        adata_bp,
        group='T_cell_fate' if color == '#E64B35' else 'Myeloid_fate',
        key='branch_deg'
    )
    # volcano: log2FC vs -log10(pval_adj)
    x = df_all['logfoldchanges'].values
    y = -np.log10(df_all['pvals_adj'].values + 1e-300)
    sig = (df_all['pvals_adj'] < 0.05) & (df_all['logfoldchanges'] > 0.5)
    ax.scatter(x[~sig], y[~sig], s=3, alpha=0.3, c='#cccccc', linewidths=0)
    ax.scatter(x[sig],  y[sig],  s=5, alpha=0.7, c=color, linewidths=0)
    # label top 10
    top10 = df_all[sig].nlargest(10, 'logfoldchanges')
    for _, row in top10.iterrows():
        ax.annotate(row['names'],
                    (row['logfoldchanges'], -np.log10(row['pvals_adj'] + 1e-300)),
                    fontsize=7, color='#222',
                    xytext=(3, 2), textcoords='offset points')
    ax.axvline(0.5, color='#888', lw=0.8, ls='--')
    ax.axhline(2,   color='#888', lw=0.8, ls='--')
    ax.set_xlabel('log2 Fold Change', fontsize=9)
    ax.set_ylabel('-log10(adj p-value)', fontsize=9)
    ax.set_title(title, fontsize=10, fontweight='bold', color=color)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('figures/phase2/branch_DEG_volcano.png',
            dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: figures/phase2/branch_DEG_volcano.png")

# ── D: UMAP showing top branch genes in branch-point zone ────────
highlight_genes = deg_T['names'].head(4).tolist() + deg_My['names'].head(4).tolist()
valid_genes = [g for g in highlight_genes if g in adata.var_names][:8]

fig, axes = plt.subplots(2, 4, figsize=(16, 8))
fig.suptitle('Top Branch-Point Genes — expression on UMAP (branch-point cells highlighted)',
             fontsize=11, fontweight='bold')
umap = adata.obsm['X_umap']
DOT_S, DOT_A = 0.5, 0.35

for ax, gene in zip(axes.flat, valid_genes):
    gene_expr = adata.obs_vector(gene)   # log-normalized after step 2
    sc_g = ax.scatter(umap[:, 0], umap[:, 1],
                      c=gene_expr, s=DOT_S, alpha=DOT_A,
                      cmap='Reds', linewidths=0, vmin=0)
    plt.colorbar(sc_g, ax=ax, shrink=0.7)
    # Overlay fate-committed cells outline
    ax.scatter(umap[keep, 0], umap[keep, 1],
               s=0.3, c='navy', alpha=0.15, linewidths=0)
    fate_side = 'T-fate' if gene in deg_T['names'].values else 'My-fate'
    ax.set_title(f'{gene}  [{fate_side}]', fontsize=9, fontweight='bold')
    ax.axis('off')

for ax in axes.flat[len(valid_genes):]:
    ax.axis('off')

plt.tight_layout()
plt.savefig('figures/phase2/branch_genes_umap.png',
            dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: figures/phase2/branch_genes_umap.png")

print(f"\n{'='*65}")
print(f"Done  |  {datetime.now().strftime('%H:%M:%S')}")
print(f"\nTop T-cell fate branch genes:")
print(deg_T[['names','logfoldchanges','pvals_adj']].head(15).to_string(index=False))
print(f"\nTop Myeloid fate branch genes:")
print(deg_My[['names','logfoldchanges','pvals_adj']].head(15).to_string(index=False))
print(f"{'='*65}")
