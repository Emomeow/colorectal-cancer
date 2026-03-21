#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 2 — Pseudotime Analysis: Pre-driver Suppressive Subpopulations
Data: GSE125527 (UC, inflammatory) + GSE132465 (CRC, tumor)
Goal: Trace immune cell evolution from inflammation → tumor microenvironment
Tool: Palantir (diffusion pseudotime on scVI latent space)
"""
import scanpy as sc
import numpy as np
import pandas as pd
import palantir
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os, sys, io
from datetime import datetime

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

os.makedirs('figures/phase2', exist_ok=True)
os.makedirs('phase2_results', exist_ok=True)

print("=" * 65)
print(f"Phase 2 — Pseudotime (Pre-driver analysis)  |  {datetime.now().strftime('%H:%M:%S')}")
print("=" * 65)

# ─────────────────────────────────────────────────────────────────
# 1. Load and subset data
# ─────────────────────────────────────────────────────────────────
print("\n[1/6] Loading and subsetting data ...")
adata_all = sc.read_h5ad('processed_data/merged_annotated.h5ad')

# Keep only GSE125527 (UC) and GSE132465 (CRC)
batch_mask = adata_all.obs['batch'].isin(['GSE125527', 'GSE132465'])

# Focus on immune cell types most relevant to immunosuppression trajectory
IMMUNE_TYPES = {'T cells', 'NK cells', 'Myeloids'}
type_mask  = adata_all.obs['cell_type_grouped'].isin(IMMUNE_TYPES)

adata = adata_all[batch_mask & type_mask].copy()
print(f"  Subset: {adata.n_obs:,} cells")
print("  Batch × Cell type breakdown:")
ct_tbl = pd.crosstab(adata.obs['batch'], adata.obs['cell_type_grouped'])
print(ct_tbl.to_string())

# ─────────────────────────────────────────────────────────────────
# 2. Build neighbourhood graph on scVI latent space
# ─────────────────────────────────────────────────────────────────
print("\n[2/6] Building KNN graph on scVI latent space ...")
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=30, random_state=0)
sc.tl.umap(adata, min_dist=0.3, random_state=0)
print(f"  UMAP computed: {adata.obsm['X_umap'].shape}")

# ─────────────────────────────────────────────────────────────────
# 3. Diffusion map (Palantir)
# ─────────────────────────────────────────────────────────────────
print("\n[3/6] Computing diffusion map ...")
dm_res = palantir.utils.run_diffusion_maps(
    pd.DataFrame(adata.obsm['X_scvi'], index=adata.obs_names),
    n_components=10,
    knn=30,
)
adata.obsm['X_diffmap'] = dm_res['EigenVectors'].values
ms_data = palantir.utils.determine_multiscale_space(dm_res)
print(f"  Diffusion map: {adata.obsm['X_diffmap'].shape}")

# ─────────────────────────────────────────────────────────────────
# 4. Choose root cell
#
#  Root = UC (GSE125527) T cell with earliest / most naive-like state.
#  Strategy: pick the cell with the smallest 1st diffusion component
#  among UC T cells — this places the root in the least differentiated
#  part of the UC compartment.
# ─────────────────────────────────────────────────────────────────
print("\n[4/6] Selecting root cell (UC T cells, naive-like end) ...")
uc_t_mask = (adata.obs['batch'] == 'GSE125527') & \
            (adata.obs['cell_type_grouped'] == 'T cells')
dc1 = dm_res['EigenVectors'].iloc[:, 0]   # 1st diffusion component
root_idx  = dc1[uc_t_mask].idxmin()       # extreme end in DC1
root_cell = root_idx
print(f"  Root cell: {root_cell}")
print(f"    batch={adata.obs.loc[root_cell,'batch']}, "
      f"type={adata.obs.loc[root_cell,'cell_type_grouped']}")

# ─────────────────────────────────────────────────────────────────
# 5. Run Palantir pseudotime
# ─────────────────────────────────────────────────────────────────
print("\n[5/6] Running Palantir pseudotime ...")
pr_res = palantir.core.run_palantir(
    ms_data,
    early_cell=root_cell,
    num_waypoints=1200,
    use_early_cell_as_start=True,
    n_jobs=4,
)
adata.obs['pseudotime']       = pr_res.pseudotime
adata.obs['entropy']          = pr_res.entropy
# branch probabilities (one column per terminal state)
for col in pr_res.branch_probs.columns:
    adata.obs[f'branch_{col}'] = pr_res.branch_probs[col].values

n_terminals = pr_res.branch_probs.shape[1]
print(f"  Pseudotime range: {pr_res.pseudotime.min():.3f} – {pr_res.pseudotime.max():.3f}")
print(f"  Terminal states identified: {n_terminals}")
print(f"  Terminal cell IDs: {list(pr_res.branch_probs.columns)}")

# Annotate terminal state identities
terminal_annots = {}
for t_cell in pr_res.branch_probs.columns:
    ct   = adata.obs.loc[t_cell, 'cell_type_grouped']
    bat  = adata.obs.loc[t_cell, 'batch']
    terminal_annots[t_cell] = f'{ct} ({bat.replace("GSE","GSE ")})'
    print(f"    {t_cell}: {terminal_annots[t_cell]}")

# ─────────────────────────────────────────────────────────────────
# 6. Plots
# ─────────────────────────────────────────────────────────────────
print("\n[6/6] Generating plots ...")

BATCH_COLORS = {'GSE125527': '#00A087', 'GSE132465': '#4DBBD5'}
TYPE_COLORS  = {'T cells': '#E64B35', 'NK cells': '#3C5488', 'Myeloids': '#00A087'}
DOT_S, DOT_A = 0.8, 0.5

umap = adata.obsm['X_umap']
rng  = np.random.default_rng(0)
shuf = rng.permutation(adata.n_obs)

# ── Fig 1: overview (batch | cell type | pseudotime | entropy) ──
fig, axes = plt.subplots(2, 2, figsize=(13, 11))
fig.suptitle('Phase 2 — Pseudotime: UC → CRC Immune Trajectory\n'
             'GSE125527 (UC) + GSE132465 (CRC)  |  T cells + NK + Myeloids',
             fontsize=12, fontweight='bold', y=1.01)

# batch
ax = axes[0, 0]
for bat, col in BATCH_COLORS.items():
    m = adata.obs['batch'].values == bat
    ax.scatter(umap[m, 0], umap[m, 1], c=col, s=DOT_S, alpha=DOT_A, linewidths=0, label=bat)
ax.set_title('By dataset (batch)', fontsize=10, fontweight='bold')
ax.axis('off')
ax.legend(markerscale=5, fontsize=8, loc='lower left', framealpha=0.7)

# cell type
ax = axes[0, 1]
for ct, col in TYPE_COLORS.items():
    m = adata.obs['cell_type_grouped'].values == ct
    if m.sum() == 0: continue
    ax.scatter(umap[m, 0], umap[m, 1], c=col, s=DOT_S, alpha=DOT_A, linewidths=0, label=ct)
ax.set_title('By cell type', fontsize=10, fontweight='bold')
ax.axis('off')
ax.legend(markerscale=5, fontsize=8, loc='lower left', framealpha=0.7)

# pseudotime
ax = axes[1, 0]
pt_vals = adata.obs['pseudotime'].values
sc_pt = ax.scatter(umap[shuf, 0], umap[shuf, 1],
                   c=pt_vals[shuf], s=DOT_S, alpha=0.6,
                   cmap='viridis', linewidths=0)
plt.colorbar(sc_pt, ax=ax, shrink=0.8, label='Pseudotime')
# Mark root
root_pos = umap[adata.obs_names.get_loc(root_cell)]
ax.scatter(*root_pos, s=120, c='red', marker='*', zorder=10, label='Root')
ax.set_title('Pseudotime (UC root → CRC terminals)', fontsize=10, fontweight='bold')
ax.axis('off')
ax.legend(fontsize=8, loc='lower left')

# entropy (branching uncertainty)
ax = axes[1, 1]
ent_vals = adata.obs['entropy'].values
sc_ent = ax.scatter(umap[shuf, 0], umap[shuf, 1],
                    c=ent_vals[shuf], s=DOT_S, alpha=0.6,
                    cmap='magma', linewidths=0)
plt.colorbar(sc_ent, ax=ax, shrink=0.8, label='Entropy')
ax.set_title('Branch entropy (high = branch point)', fontsize=10, fontweight='bold')
ax.axis('off')

plt.tight_layout()
out1 = 'figures/phase2/pseudotime_overview.png'
plt.savefig(out1, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: {out1}")

# ── Fig 2: branch probabilities (one panel per terminal) ──────────
fig, axes = plt.subplots(1, n_terminals, figsize=(5 * n_terminals, 5))
if n_terminals == 1:
    axes = [axes]
fig.suptitle('Branch Probabilities toward Terminal States',
             fontsize=11, fontweight='bold', y=1.01)

for ax, t_cell in zip(axes, pr_res.branch_probs.columns):
    bp = pr_res.branch_probs[t_cell].values
    sc_bp = ax.scatter(umap[shuf, 0], umap[shuf, 1],
                       c=bp[shuf], s=DOT_S, alpha=0.6,
                       cmap='YlOrRd', vmin=0, vmax=1, linewidths=0)
    plt.colorbar(sc_bp, ax=ax, shrink=0.8)
    ax.set_title(f'→ {terminal_annots[t_cell]}', fontsize=9, fontweight='bold')
    ax.axis('off')

plt.tight_layout()
out2 = 'figures/phase2/pseudotime_branch_probs.png'
plt.savefig(out2, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: {out2}")

# ── Fig 3: pseudotime distribution per batch × cell type ──────────
fig, ax = plt.subplots(figsize=(9, 4))
groups = {
    'UC T cells'      : (adata.obs['batch']=='GSE125527') & (adata.obs['cell_type_grouped']=='T cells'),
    'UC NK'           : (adata.obs['batch']=='GSE125527') & (adata.obs['cell_type_grouped']=='NK cells'),
    'UC Myeloids'     : (adata.obs['batch']=='GSE125527') & (adata.obs['cell_type_grouped']=='Myeloids'),
    'CRC T cells'     : (adata.obs['batch']=='GSE132465') & (adata.obs['cell_type_grouped']=='T cells'),
    'CRC NK'          : (adata.obs['batch']=='GSE132465') & (adata.obs['cell_type_grouped']=='NK cells'),
    'CRC Myeloids'    : (adata.obs['batch']=='GSE132465') & (adata.obs['cell_type_grouped']=='Myeloids'),
}
colors_vio = ['#00A087','#3C5488','#F39B7F','#E64B35','#4DBBD5','#B09C85']
data_vio   = [adata.obs['pseudotime'].values[m.values] for m in groups.values()]
vp = ax.violinplot(data_vio, positions=range(len(groups)), showmedians=True)
for i, (body, col) in enumerate(zip(vp['bodies'], colors_vio)):
    body.set_facecolor(col); body.set_alpha(0.7)
ax.set_xticks(range(len(groups)))
ax.set_xticklabels(list(groups.keys()), rotation=30, ha='right', fontsize=9)
ax.set_ylabel('Pseudotime', fontsize=10)
ax.set_title('Pseudotime distribution: UC vs CRC across cell types', fontsize=11, fontweight='bold')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
plt.tight_layout()
out3 = 'figures/phase2/pseudotime_distribution.png'
plt.savefig(out3, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  Saved: {out3}")

# ── Save results ──────────────────────────────────────────────────
pr_res.pseudotime.to_csv('phase2_results/pseudotime_values.csv')
pr_res.branch_probs.to_csv('phase2_results/branch_probs.csv')
adata.obs[['batch','cell_type_grouped','pseudotime','entropy']].to_csv(
    'phase2_results/cell_metadata_with_pseudotime.csv')
print("  Saved: phase2_results/pseudotime_values.csv")
print("  Saved: phase2_results/branch_probs.csv")
print("  Saved: phase2_results/cell_metadata_with_pseudotime.csv")

# Save adata for downstream gene trend analysis
adata.write_h5ad('phase2_results/adata_pseudotime.h5ad')
print("  Saved: phase2_results/adata_pseudotime.h5ad")

print(f"\n{'='*65}")
print(f"Done  |  {datetime.now().strftime('%H:%M:%S')}")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Terminal states: {n_terminals}")
print(f"  Outputs: figures/phase2/  +  phase2_results/")
print(f"{'='*65}")
