#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Replot gene trend heatmap at reasonable size"""
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, io, os

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

os.makedirs('figures/phase2', exist_ok=True)

print("Loading adata_hvg from saved pseudotime results ...")
import palantir

# Reload HVG adata and recompute gene trends (fast, ~5s)
adata = sc.read_h5ad('phase2_results/adata_pseudotime.h5ad')
pr_pseudotime   = pd.read_csv('phase2_results/pseudotime_values.csv', index_col=0).squeeze()
pr_branch_probs = pd.read_csv('phase2_results/branch_probs.csv', index_col=0)

adata.obs['palantir_pseudotime'] = pr_pseudotime.reindex(adata.obs_names).values
adata.obsm['palantir_fate_probabilities'] = pr_branch_probs.reindex(adata.obs_names)

T_TERM  = pr_branch_probs.columns[0]
MY_TERM = pr_branch_probs.columns[1]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2500, batch_key='batch')
adata_hvg = adata[:, adata.var['highly_variable']].copy()
adata_hvg.obs['palantir_pseudotime'] = adata.obs['palantir_pseudotime']
adata_hvg.obsm['palantir_fate_probabilities'] = adata.obsm['palantir_fate_probabilities']

print("Computing gene trends ...")
palantir.presults.select_branch_cells(adata_hvg, q=0.01, eps=0.01)
gene_trends = palantir.presults.compute_gene_trends(
    adata_hvg, lineages=[T_TERM, MY_TERM], pseudo_time_key='palantir_pseudotime',
)

# ── Extract trend matrices ────────────────────────────────────────
# compute_gene_trends returns a dict keyed by terminal cell IDs
print(f"Gene trends keys: {list(gene_trends.keys())[:2]}")

df_T  = gene_trends[T_TERM]['trends']
df_My = gene_trends[MY_TERM]['trends']

if df_T is None or df_My is None:
    print("ERROR: trend data not found")
    sys.exit(1)

print(f"T branch trends:  {df_T.shape}  (genes × pseudotime_points)")
print(f"My branch trends: {df_My.shape}")

# ── Load DEG results to get top genes ────────────────────────────
deg_T  = pd.read_csv('phase2_results/branch_DEG_Tcell_fate.csv')
deg_My = pd.read_csv('phase2_results/branch_DEG_Myeloid_fate.csv')

# Gene names are in the index (rows) of the trends DataFrames
trend_gene_set = set(df_T.index) & set(df_My.index)
top_T_genes  = [g for g in deg_T['names'].head(50)  if g in trend_gene_set][:25]
top_My_genes = [g for g in deg_My['names'].head(50) if g in trend_gene_set][:25]

# Fallback: use all trend genes ranked by variance if DEG overlap is low
if len(top_T_genes) < 10:
    var_T = df_T.var(axis=1).sort_values(ascending=False)
    top_T_genes = var_T.head(25).index.tolist()
if len(top_My_genes) < 10:
    var_My = df_My.var(axis=1).sort_values(ascending=False)
    top_My_genes = var_My.head(25).index.tolist()

print(f"Top T genes for heatmap:  {len(top_T_genes)}")
print(f"Top My genes for heatmap: {len(top_My_genes)}")

# ── Build heatmap matrices ────────────────────────────────────────
# Downsample pseudotime axis to 100 points
N_PT = 100
def resample_trends(df, genes, n_pts=N_PT):
    """Resample trend DataFrame to n_pts pseudotime bins, z-score normalize."""
    df_sel = df.loc[genes]   # genes × pseudotime_points (rows=genes, cols=pseudotime)
    # Resample to n_pts columns
    idx = np.linspace(0, df_sel.shape[1] - 1, n_pts).astype(int)
    df_rs = df_sel.iloc[:, idx]
    # Z-score per gene
    mu   = df_rs.mean(axis=1)
    sig  = df_rs.std(axis=1).clip(lower=1e-6)
    df_z = df_rs.subtract(mu, axis=0).divide(sig, axis=0)
    return df_z.clip(-3, 3)

mat_T_T   = resample_trends(df_T,  top_T_genes)   # T-branch genes on T lineage
mat_T_My  = resample_trends(df_T,  top_My_genes)  # My-branch genes on T lineage
mat_My_T  = resample_trends(df_My, top_T_genes)   # T-branch genes on My lineage
mat_My_My = resample_trends(df_My, top_My_genes)  # My-branch genes on My lineage

# ── Plot: 2×2 heatmap grid ────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 10),
                          gridspec_kw={'hspace': 0.35, 'wspace': 0.4})
fig.suptitle('Gene Expression Trends Along Pseudotime Lineages\n'
             'T-cell branch (left col) vs Myeloid branch (right col)',
             fontsize=12, fontweight='bold')

CMAP = 'RdBu_r'
VMAX = 2.5

panels = [
    (axes[0, 0], mat_T_T,   top_T_genes,  'T-branch genes  |  T-cell lineage',  '#E64B35'),
    (axes[0, 1], mat_My_My, top_My_genes, 'My-branch genes  |  Myeloid lineage', '#00A087'),
    (axes[1, 0], mat_T_My,  top_My_genes, 'My-branch genes  |  T-cell lineage',  '#4DBBD5'),
    (axes[1, 1], mat_My_T,  top_T_genes,  'T-branch genes  |  Myeloid lineage',  '#F39B7F'),
]

for ax, mat, genes, title, title_col in panels:
    im = ax.imshow(mat.values, aspect='auto', cmap=CMAP,
                   vmin=-VMAX, vmax=VMAX, interpolation='nearest')
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=7)
    ax.set_xticks([0, N_PT // 2, N_PT - 1])
    ax.set_xticklabels(['Early\n(UC)', 'Mid', 'Late\n(CRC)'], fontsize=8)
    ax.set_title(title, fontsize=9, fontweight='bold', color=title_col, pad=4)
    plt.colorbar(im, ax=ax, shrink=0.8, label='Z-score', pad=0.01)

plt.savefig('figures/phase2/gene_trend_heatmaps.png',
            dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print("Saved: figures/phase2/gene_trend_heatmaps.png")

# ── Also save trend data as CSV for further analysis ──────────────
df_T.loc[top_T_genes + top_My_genes].to_csv('phase2_results/gene_trends_T_lineage.csv')
df_My.loc[top_T_genes + top_My_genes].to_csv('phase2_results/gene_trends_My_lineage.csv')
print("Saved: phase2_results/gene_trends_*.csv")
print("Done.")
