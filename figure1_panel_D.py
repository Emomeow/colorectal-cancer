#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 1 Panel D — Marker gene dotplot across all cell types
"""
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os, sys, io

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

os.makedirs('figures', exist_ok=True)

print("Loading data ...")
adata = sc.read_h5ad('processed_data/merged_annotated.h5ad')
print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

# ── Normalise for visualisation (work on a copy) ──────────────────────────
adata_vis = adata.copy()
sc.pp.normalize_total(adata_vis, target_sum=1e4)
sc.pp.log1p(adata_vis)

# ── Marker gene dictionary ────────────────────────────────────────────────
# Ordered by biological lineage; only genes confirmed present in panel
marker_dict = {
    # ── Immune: T cells ──────────────────────────────────────────────
    'T cells\n(pan)'       : ['CD3D', 'CD3E', 'CD2'],
    'CD4 T'                : ['CD4', 'IL7R', 'CCR7', 'TCF7'],
    'CD8 T'                : ['CD8A', 'CD8B'],
    'Treg'                 : ['FOXP3', 'TIGIT', 'IL2RA'],
    'Effector / Exhaust'   : ['GZMB', 'GZMK', 'PDCD1', 'LAG3', 'HAVCR2'],
    # ── Immune: NK cells ─────────────────────────────────────────────
    'NK cells'             : ['NCAM1', 'KLRD1', 'KLRB1', 'KLRC1', 'XCL1'],
    # ── Immune: B cells ──────────────────────────────────────────────
    'B cells\n(pan)'       : ['CD79A', 'MS4A1', 'CD19'],
    'Plasma B'             : ['CD38', 'IGHG1', 'IGHA1', 'IGKC', 'MZB1'],
    # ── Immune: Myeloids ─────────────────────────────────────────────
    'Monocyte'             : ['CD14', 'LYZ', 'S100A8', 'S100A9'],
    'Macrophage'           : ['CD68', 'CSF1R', 'CD163', 'FCGR3A'],
    # ── Non-immune ───────────────────────────────────────────────────
    'Epithelial'           : ['KRT18'],
    'Malignant'            : ['ASCL2', 'PCNA'],
    'Stromal'              : ['ACTA2', 'VIM'],
    'Mast cells'           : ['FCER1A'],
}

# Flatten and verify all genes present
all_markers = []
marker_dict_clean = {}
for group, genes in marker_dict.items():
    available = [g for g in genes if g in adata_vis.var_names]
    if available:
        marker_dict_clean[group] = available
        all_markers.extend(available)
all_markers = list(dict.fromkeys(all_markers))  # deduplicate, preserve order

print(f"  Marker groups: {len(marker_dict_clean)}")
print(f"  Total marker genes: {len(all_markers)}")
for grp, genes in marker_dict_clean.items():
    print(f"    {grp:22s}: {genes}")

# ── Cell type order for y-axis ────────────────────────────────────────────
ct_order = [
    'T cells',
    'NK cells',
    'B cells',
    'Myeloids',
    'Stromal cells',
    'Mast cells',
    'Epithelial cells',
    'BEST4+ colonocyte',
    'colon goblet cell',
    'colonocyte',
    'early colonocyte',
    'secretory cell',
    'intestinal crypt stem cell of colon',
    'enteroendocrine cell of colon',
    'tuft cell of colon',
    'malignant cell',
    'Unknown',
]
# Keep only types present in data
ct_present = [ct for ct in ct_order
              if ct in adata_vis.obs['unified_cell_type'].cat.categories]
# Add any remaining types not in our order list
extra = [ct for ct in adata_vis.obs['unified_cell_type'].cat.categories
         if ct not in ct_present and ct not in ['Unspecified']]
ct_present = ct_present + extra

print(f"\nCell types in plot ({len(ct_present)}): {ct_present}")

# Filter to non-trivial types (>10 cells)
ct_present = [ct for ct in ct_present
              if (adata_vis.obs['unified_cell_type'] == ct).sum() > 10]

# ── Compute dot values manually for full control ──────────────────────────
print("\nComputing dot values ...")

def compute_dots(adata, cell_types, genes):
    """Returns mean_expr and frac_expressing DataFrames."""
    means = pd.DataFrame(index=cell_types, columns=genes, dtype=float)
    fracs = pd.DataFrame(index=cell_types, columns=genes, dtype=float)
    for ct in cell_types:
        mask = adata.obs['unified_cell_type'] == ct
        X_ct = adata[mask, genes].X
        if hasattr(X_ct, 'toarray'):
            X_ct = X_ct.toarray()
        means.loc[ct] = X_ct.mean(axis=0)
        fracs.loc[ct] = (X_ct > 0).mean(axis=0)
    return means, fracs

means, fracs = compute_dots(adata_vis, ct_present, all_markers)

# Scale means per gene to [0,1] for colour (so each gene uses full range)
means_scaled = means.copy()
for gene in all_markers:
    col = means[gene]
    rng = col.max() - col.min()
    means_scaled[gene] = (col - col.min()) / rng if rng > 0 else 0

# ── Plot ──────────────────────────────────────────────────────────────────
print("Plotting ...")

n_ct   = len(ct_present)
n_gene = len(all_markers)

fig_w  = max(16, n_gene * 0.42 + 3)
fig_h  = max(8,  n_ct   * 0.42 + 3)
# 右侧留出固定宽度给 colorbar + size legend
RMARGIN = 1.8   # inches
fig    = plt.figure(figsize=(fig_w + RMARGIN, fig_h))

# 主绘图区：左10%，右侧留给图例
ax_left   = 0.08
ax_bottom = 0.16
ax_width  = (fig_w - 0.5) / (fig_w + RMARGIN)
ax_height = 0.74
ax = fig.add_axes([ax_left, ax_bottom, ax_width, ax_height])

# 计算右侧两个图例的 figure 坐标
r0 = ax_left + ax_width + 0.02          # 右侧起点 x（figure 坐标）
cbar_x  = r0 + 0.01
cbar_y  = ax_bottom + ax_height * 0.52
cbar_w  = 0.018
cbar_h  = ax_height * 0.38
leg_x   = r0
leg_y   = ax_bottom + ax_height * 0.05
leg_w   = RMARGIN / (fig_w + RMARGIN) - 0.03
leg_h   = ax_height * 0.42

# Draw dots
cmap = matplotlib.colormaps['Reds']
max_dot = 10   # pt radius scaling

for i, ct in enumerate(ct_present):
    for j, gene in enumerate(all_markers):
        color   = cmap(float(means_scaled.loc[ct, gene]))
        frac    = float(fracs.loc[ct, gene])
        size    = (frac ** 0.5) * max_dot   # area ∝ fraction
        ax.scatter(j, i, s=size**2, c=[color],
                   linewidths=0.3, edgecolors='#888888', zorder=2)

# Grid lines
ax.set_xticks(range(n_gene))
ax.set_xticklabels(all_markers, rotation=90, fontsize=8)
ax.set_yticks(range(n_ct))
ax.set_yticklabels(ct_present, fontsize=9)
ax.set_xlim(-0.7, n_gene - 0.3)
ax.set_ylim(-0.7, n_ct  - 0.3)
ax.grid(True, linewidth=0.3, alpha=0.4, zorder=1)
ax.set_axisbelow(True)

# Marker group separators + labels on top
gene_pos = 0
for grp, genes in marker_dict_clean.items():
    n = len(genes)
    mid = gene_pos + n / 2 - 0.5
    ax.text(mid, 1.02, grp, ha='center', va='bottom',
            fontsize=7.5, fontweight='bold', color='#333333',
            transform=ax.get_xaxis_transform())
    if gene_pos > 0:
        ax.axvline(gene_pos - 0.5, color='#aaaaaa', linewidth=0.8, zorder=3)
    gene_pos += n

# ── Colorbar — fig.add_axes 精确定位，绝不进入主绘图区 ────────────────────
from matplotlib.lines import Line2D

cax = fig.add_axes([cbar_x, cbar_y, cbar_w, cbar_h])
sm  = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
sm.set_array([])
cbar = fig.colorbar(sm, cax=cax)
cbar.set_label('Scaled mean\nexpression', fontsize=8)
cbar.ax.tick_params(labelsize=7)

# ── Dot size legend — 专用 axes，手动绘制 ────────────────────────────────
leg_ax = fig.add_axes([leg_x, leg_y, leg_w, leg_h])
leg_ax.set_xlim(0, 1); leg_ax.set_ylim(0, 1); leg_ax.axis('off')
fracs_leg = [1.0, 0.75, 0.5, 0.25, 0.1]
n_leg = len(fracs_leg)
for k, frac in enumerate(fracs_leg):
    y = 0.82 - k * (0.72 / (n_leg - 1))
    s = (frac ** 0.5) * max_dot
    leg_ax.scatter(0.18, y, s=s**2, c='grey',
                   linewidths=0.3, edgecolors='#555555',
                   transform=leg_ax.transData, zorder=5)
    leg_ax.text(0.35, y, f'{int(frac*100)}%', va='center',
                fontsize=7.5, transform=leg_ax.transData)
leg_ax.text(0.18, 0.96, '% expressing', va='top', ha='center',
            fontsize=8, fontweight='bold', transform=leg_ax.transData)
leg_ax.patch.set_facecolor('white')
leg_ax.patch.set_alpha(0.7)
for spine in leg_ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.5)
    spine.set_edgecolor('#aaaaaa')

ax.set_title('Figure 1D — Marker Gene Expression across Cell Types',
             fontsize=13, fontweight='bold', pad=28)
ax.set_xlabel('Marker genes', fontsize=10, labelpad=6)
ax.set_ylabel('Cell type', fontsize=10)

out = 'figures/figure1_panel_D_dotplot.png'
plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved: {out}")

# Also save as PDF for publication
out_pdf = 'figures/figure1_panel_D_dotplot.pdf'
fig2 = plt.figure(figsize=(fig_w + RMARGIN, fig_h))
ax2  = fig2.add_axes([ax_left, ax_bottom, ax_width, ax_height])

for i, ct in enumerate(ct_present):
    for j, gene in enumerate(all_markers):
        color = cmap(float(means_scaled.loc[ct, gene]))
        frac  = float(fracs.loc[ct, gene])
        size  = (frac ** 0.5) * max_dot
        ax2.scatter(j, i, s=size**2, c=[color],
                    linewidths=0.3, edgecolors='#888888', zorder=2)
ax2.set_xticks(range(n_gene))
ax2.set_xticklabels(all_markers, rotation=90, fontsize=8)
ax2.set_yticks(range(n_ct))
ax2.set_yticklabels(ct_present, fontsize=9)
ax2.set_xlim(-0.7, n_gene - 0.3)
ax2.set_ylim(-0.7, n_ct  - 0.3)
ax2.grid(True, linewidth=0.3, alpha=0.4, zorder=1)
ax2.set_axisbelow(True)
gene_pos = 0
for grp, genes in marker_dict_clean.items():
    n = len(genes)
    mid = gene_pos + n / 2 - 0.5
    ax2.text(mid, 1.02, grp, ha='center', va='bottom',
             fontsize=7.5, fontweight='bold', color='#333333',
             transform=ax2.get_xaxis_transform())
    if gene_pos > 0:
        ax2.axvline(gene_pos - 0.5, color='#aaaaaa', linewidth=0.8, zorder=3)
    gene_pos += n

# Colorbar — fig.add_axes for precise placement outside main axes
cax2 = fig2.add_axes([cbar_x, cbar_y, cbar_w, cbar_h])
sm2  = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
sm2.set_array([])
cbar2 = fig2.colorbar(sm2, cax=cax2)
cbar2.set_label('Scaled mean\nexpression', fontsize=8)
cbar2.ax.tick_params(labelsize=7)

# Dot size legend — dedicated axes, manual scatter
leg_ax2 = fig2.add_axes([leg_x, leg_y, leg_w, leg_h])
leg_ax2.set_xlim(0, 1); leg_ax2.set_ylim(0, 1); leg_ax2.axis('off')
for k, frac in enumerate(fracs_leg):
    y = 0.82 - k * (0.72 / (n_leg - 1))
    s = (frac ** 0.5) * max_dot
    leg_ax2.scatter(0.18, y, s=s**2, c='grey',
                    linewidths=0.3, edgecolors='#555555',
                    transform=leg_ax2.transData, zorder=5)
    leg_ax2.text(0.35, y, f'{int(frac*100)}%', va='center',
                 fontsize=7.5, transform=leg_ax2.transData)
leg_ax2.text(0.18, 0.96, '% expressing', va='top', ha='center',
             fontsize=8, fontweight='bold', transform=leg_ax2.transData)
leg_ax2.patch.set_facecolor('white')
leg_ax2.patch.set_alpha(0.7)
for spine in leg_ax2.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.5)
    spine.set_edgecolor('#aaaaaa')

ax2.set_title('Figure 1D — Marker Gene Expression across Cell Types',
              fontsize=13, fontweight='bold', pad=28)
ax2.set_xlabel('Marker genes', fontsize=10, labelpad=6)
ax2.set_ylabel('Cell type', fontsize=10)

plt.savefig(out_pdf, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved: {out_pdf}")

# Copy both outputs to phase1 folder
import shutil
os.makedirs('figures/phase1', exist_ok=True)
shutil.copy(out,     'figures/phase1/figure1_panel_D_dotplot.png')
shutil.copy(out_pdf, 'figures/phase1/figure1_panel_D_dotplot.pdf')
print("Copied to figures/phase1/")
print("Done.")
