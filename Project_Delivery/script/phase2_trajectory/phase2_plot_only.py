#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Re-plot pseudotime results from saved h5ad (skip recomputation)"""
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys, io

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

os.makedirs('figures/phase2', exist_ok=True)

print("Loading saved pseudotime data ...")
adata = sc.read_h5ad('phase2_results/adata_pseudotime.h5ad')
pr_pseudotime  = pd.read_csv('phase2_results/pseudotime_values.csv', index_col=0).squeeze()
pr_branch_probs = pd.read_csv('phase2_results/branch_probs.csv', index_col=0)

umap = adata.obsm['X_umap']
rng  = np.random.default_rng(0)
shuf = rng.permutation(adata.n_obs)

BATCH_COLORS = {'GSE125527': '#00A087', 'GSE132465': '#4DBBD5'}
TYPE_COLORS  = {'T cells': '#E64B35', 'NK cells': '#3C5488', 'Myeloids': '#00A087'}
DOT_S, DOT_A = 0.8, 0.45

# terminal annotations
terminal_annots = {}
for t_cell in pr_branch_probs.columns:
    ct  = adata.obs.loc[t_cell, 'cell_type_grouped']
    bat = adata.obs.loc[t_cell, 'batch']
    terminal_annots[t_cell] = f'{ct}\n({bat})'
n_terminals = len(pr_branch_probs.columns)
print(f"  {adata.n_obs:,} cells  |  {n_terminals} terminal states")
for k, v in terminal_annots.items():
    print(f"    {v.replace(chr(10),' ')}")

# ── Fig 1: 4-panel overview ────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(13, 11))
fig.suptitle('Phase 2 — Pseudotime: UC (GSE125527) → CRC (GSE132465)\n'
             'T cells + NK cells + Myeloids  |  Palantir on scVI latent space',
             fontsize=12, fontweight='bold', y=1.01)

# batch
ax = axes[0, 0]
for bat, col in BATCH_COLORS.items():
    m = adata.obs['batch'].values == bat
    ax.scatter(umap[shuf[m[shuf]], 0], umap[shuf[m[shuf]], 1],
               c=col, s=DOT_S, alpha=DOT_A, linewidths=0, label=bat)
ax.set_title('By dataset', fontsize=10, fontweight='bold')
ax.axis('off')
ax.legend(markerscale=5, fontsize=8, loc='lower left', framealpha=0.8)

# cell type
ax = axes[0, 1]
for ct, col in TYPE_COLORS.items():
    m = adata.obs['cell_type_grouped'].values == ct
    if m.sum() == 0: continue
    ax.scatter(umap[m, 0], umap[m, 1], c=col, s=DOT_S, alpha=DOT_A, linewidths=0, label=ct)
ax.set_title('By cell type', fontsize=10, fontweight='bold')
ax.axis('off')
ax.legend(markerscale=5, fontsize=8, loc='lower left', framealpha=0.8)

# pseudotime
ax = axes[1, 0]
pt = pr_pseudotime.reindex(adata.obs_names).values
sc_pt = ax.scatter(umap[shuf, 0], umap[shuf, 1],
                   c=pt[shuf], s=DOT_S, alpha=0.6,
                   cmap='viridis', linewidths=0)
plt.colorbar(sc_pt, ax=ax, shrink=0.8, label='Pseudotime')
root_mask = pr_pseudotime.reindex(adata.obs_names) == pr_pseudotime.min()
root_idx  = np.where(root_mask.values)[0][0]
ax.scatter(umap[root_idx, 0], umap[root_idx, 1],
           s=150, c='red', marker='*', zorder=10, label='Root (UC)')
for t_cell in pr_branch_probs.columns:
    tidx = adata.obs_names.get_loc(t_cell)
    ax.scatter(umap[tidx, 0], umap[tidx, 1],
               s=120, c='orange', marker='D', zorder=10)
ax.set_title('Pseudotime  (★ root · ◆ terminals)', fontsize=10, fontweight='bold')
ax.axis('off')
ax.legend(fontsize=8, loc='lower left')

# entropy
ax = axes[1, 1]
ent = adata.obs['entropy'].values
sc_ent = ax.scatter(umap[shuf, 0], umap[shuf, 1],
                    c=ent[shuf], s=DOT_S, alpha=0.6,
                    cmap='magma', linewidths=0)
plt.colorbar(sc_ent, ax=ax, shrink=0.8, label='Branch entropy')
ax.set_title('Branch entropy  (high = branching zone)', fontsize=10, fontweight='bold')
ax.axis('off')

plt.tight_layout()
out1 = 'figures/phase2/pseudotime_overview.png'
plt.savefig(out1, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved: {out1}")

# ── Fig 2: branch probabilities ───────────────────────────────────
fig, axes = plt.subplots(1, n_terminals, figsize=(5.5 * n_terminals, 5))
if n_terminals == 1: axes = [axes]
fig.suptitle('Branch Probabilities toward Terminal States',
             fontsize=11, fontweight='bold')
for ax, t_cell in zip(axes, pr_branch_probs.columns):
    bp = pr_branch_probs[t_cell].reindex(adata.obs_names).values
    sc_bp = ax.scatter(umap[shuf, 0], umap[shuf, 1],
                       c=bp[shuf], s=DOT_S, alpha=0.6,
                       cmap='YlOrRd', vmin=0, vmax=1, linewidths=0)
    plt.colorbar(sc_bp, ax=ax, shrink=0.8)
    tidx = adata.obs_names.get_loc(t_cell)
    ax.scatter(umap[tidx, 0], umap[tidx, 1], s=150, c='blue', marker='D', zorder=10)
    ax.set_title(f'→ {terminal_annots[t_cell]}', fontsize=10, fontweight='bold')
    ax.axis('off')
plt.tight_layout()
out2 = 'figures/phase2/pseudotime_branch_probs.png'
plt.savefig(out2, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved: {out2}")

# ── Fig 3: pseudotime violin per group ───────────────────────────
groups = {
    'UC T cells'  : (adata.obs['batch']=='GSE125527') & (adata.obs['cell_type_grouped']=='T cells'),
    'UC NK'       : (adata.obs['batch']=='GSE125527') & (adata.obs['cell_type_grouped']=='NK cells'),
    'UC Myeloids' : (adata.obs['batch']=='GSE125527') & (adata.obs['cell_type_grouped']=='Myeloids'),
    'CRC T cells' : (adata.obs['batch']=='GSE132465') & (adata.obs['cell_type_grouped']=='T cells'),
    'CRC NK'      : (adata.obs['batch']=='GSE132465') & (adata.obs['cell_type_grouped']=='NK cells'),
    'CRC Myeloids': (adata.obs['batch']=='GSE132465') & (adata.obs['cell_type_grouped']=='Myeloids'),
}
colors_vio = ['#00A087','#3C5488','#F39B7F','#E64B35','#4DBBD5','#B09C85']
data_vio   = [pt[m.values] for m in groups.values()]

fig, ax = plt.subplots(figsize=(10, 4.5))
vp = ax.violinplot(data_vio, positions=range(len(groups)), showmedians=True, showextrema=False)
for i, (body, col) in enumerate(zip(vp['bodies'], colors_vio)):
    body.set_facecolor(col); body.set_alpha(0.75)
vp['cmedians'].set_color('black'); vp['cmedians'].set_linewidth(1.5)
ax.set_xticks(range(len(groups)))
ax.set_xticklabels(list(groups.keys()), rotation=28, ha='right', fontsize=9)
ax.set_ylabel('Pseudotime', fontsize=10)
ax.set_title('Pseudotime distribution: UC (early) vs CRC (late) across cell types',
             fontsize=11, fontweight='bold')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
# UC/CRC brackets
ax.axvline(2.5, color='#cccccc', lw=1, ls='--')
ax.text(1, ax.get_ylim()[1]*0.97, 'UC (GSE125527)', ha='center', fontsize=8.5,
        color='#00A087', fontweight='bold')
ax.text(4, ax.get_ylim()[1]*0.97, 'CRC (GSE132465)', ha='center', fontsize=8.5,
        color='#4DBBD5', fontweight='bold')
plt.tight_layout()
out3 = 'figures/phase2/pseudotime_distribution.png'
plt.savefig(out3, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved: {out3}")
print("Done.")
