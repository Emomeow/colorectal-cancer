#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 1 – Panel D: Canonical marker dotplot (CellAssign cell types, all batches)"""

import scanpy as sc
import numpy as np
import scipy.sparse as sp
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
plt.rcParams['font.family'] = 'Arial'
os.makedirs('figures/phase1', exist_ok=True)

adata = sc.read_h5ad('scvi_results_hvg3k/merged_scvi_hvg3k.h5ad')

# ── Cell type order (CellAssign labels) ───────────────────────────────────────
CT_ORDER = ['CD4 T', 'CD8 T', 'Treg', 'NK',
            'B cell', 'Plasma',
            'Myeloid', 'Mast',
            'Stromal', 'Epithelial']
# Filter to types with > 50 cells
ct_counts = adata.obs['cellass_type'].value_counts()
CT_ORDER = [c for c in CT_ORDER if ct_counts.get(c, 0) > 50]

# ── Marker genes (confirmed present in HVG3k) ─────────────────────────────────
avail = set(adata.var_names)
MARKER_GROUPS = {
    'T / NK':   [g for g in ['CD3D', 'CD3E', 'CD4', 'CD8A', 'NCAM1', 'GNLY', 'NKG7'] if g in avail],
    'Treg':     [g for g in ['FOXP3', 'IL2RA', 'CTLA4'] if g in avail],
    'B / Plasma':[g for g in ['CD19', 'CD79A', 'MS4A1', 'MZB1', 'IGKC'] if g in avail],
    'Myeloid':  [g for g in ['LYZ', 'CD14', 'CD68', 'CSF1R'] if g in avail],
    'Epithelial':[g for g in ['KRT18', 'EPCAM', 'MUC2'] if g in avail],
    'Stroma':   [g for g in ['VIM', 'ACTA2'] if g in avail],
}
MARKER_GROUPS = {k: v for k, v in MARKER_GROUPS.items() if v}

gene_list   = []
group_spans = {}
pos = 0
for grp, genes in MARKER_GROUPS.items():
    gene_list.extend(genes)
    group_spans[grp] = (pos, pos + len(genes))
    pos += len(genes)

print(f"Cell types: {CT_ORDER}")
print(f"Genes ({len(gene_list)}): {gene_list}")

# ── Compute mean log-expression & fraction expressing ─────────────────────────
gene_idx = np.array([np.where(adata.var_names == g)[0][0] for g in gene_list])
X_raw = adata.layers['counts']
if sp.issparse(X_raw):
    X_raw = X_raw.toarray()
X_log = np.log1p(X_raw)

n_ct = len(CT_ORDER)
n_g  = len(gene_list)
dot_mean = np.zeros((n_ct, n_g))
dot_frac = np.zeros((n_ct, n_g))

for i, ct in enumerate(CT_ORDER):
    mask = (adata.obs['cellass_type'] == ct).values
    if mask.sum() < 2:
        continue
    sub_log = X_log[mask][:, gene_idx]
    sub_raw = X_raw[mask][:, gene_idx]
    dot_mean[i] = sub_log.mean(axis=0)
    dot_frac[i] = (sub_raw > 0).mean(axis=0)

# Normalize mean per gene to [0, 1]
gene_max = dot_mean.max(axis=0, keepdims=True)
gene_max[gene_max == 0] = 1.0
dot_norm = dot_mean / gene_max

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 5.0))
fig.subplots_adjust(left=0.16, right=0.97, top=0.86, bottom=0.18)

DOT_SCALE = 220
CMAP = cm.YlOrRd

for i, ct in enumerate(CT_ORDER):
    y = n_ct - 1 - i
    for j in range(n_g):
        frac = dot_frac[i, j]
        if frac < 0.02:
            continue
        cval  = dot_norm[i, j]
        color = CMAP(0.15 + cval * 0.85)
        ax.scatter(j, y, s=frac * DOT_SCALE, c=[color],
                   linewidths=0.3, edgecolors='#aaaaaa', alpha=0.95, zorder=3)

# Gene group shading + header
for grp, (x0, x1) in group_spans.items():
    if (x1 - x0) > 0:
        ax.axvspan(x0 - 0.5, x1 - 0.5, alpha=0.06, color='#888', lw=0)
        ax.text((x0 + x1) / 2 - 0.5, n_ct - 0.25, grp,
                fontsize=7, ha='center', va='bottom', color='#333', style='italic')

# Separator: immune vs stromal/epithelial
n_immune = sum(1 for c in CT_ORDER
               if c in ('CD4 T', 'CD8 T', 'Treg', 'NK', 'B cell', 'Plasma', 'Myeloid', 'Mast'))
sep_y = n_ct - n_immune - 0.5
ax.axhline(sep_y, color='#aaa', lw=0.7, ls='--', alpha=0.6)

ax.set_xlim(-0.5, n_g - 0.5)
ax.set_ylim(-0.5, n_ct + 0.2)
ax.set_xticks(range(n_g))
ax.set_xticklabels(gene_list, rotation=90, fontsize=7.5, style='italic')
ax.set_yticks(range(n_ct))
ax.set_yticklabels([CT_ORDER[n_ct - 1 - i] for i in range(n_ct)], fontsize=7.5)
ax.tick_params(labelsize=7.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.text(-0.14, 1.04, 'D', transform=ax.transAxes,
        fontsize=14, fontweight='bold', va='bottom')
ax.set_title('Canonical Marker Genes (CellAssign)', fontsize=10, pad=20)

# Color bar
import matplotlib.colors as mcolors
sm = plt.cm.ScalarMappable(cmap=CMAP, norm=mcolors.Normalize(vmin=0, vmax=1))
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, shrink=0.35, aspect=12,
                    location='right', pad=0.01)
cbar.set_label('Mean expression\n(normalized)', fontsize=7)
cbar.ax.tick_params(labelsize=6.5)

# Dot size legend — snapped to colorbar bottom
for fv, lab in [(0.1, '10%'), (0.5, '50%'), (1.0, '100%')]:
    ax.scatter([], [], s=fv * DOT_SCALE, c='#888', alpha=0.7,
               linewidths=0.3, edgecolors='#888', label=lab)
fig.canvas.draw()
cb_pos = cbar.ax.get_position()
ax_pos = ax.get_position()
cbar_bottom_in_ax = (cb_pos.y0 - ax_pos.y0) / ax_pos.height
size_leg = ax.legend(fontsize=7, frameon=False, title='% expressing',
                     title_fontsize=7, loc='upper left',
                     bbox_to_anchor=(1.01, cbar_bottom_in_ax),
                     handletextpad=0.3, labelspacing=1.2)
ax.add_artist(size_leg)

for ext in ['png', 'pdf']:
    fig.savefig(f'figures/phase1/fig1_panelD.{ext}',
                dpi=200, bbox_inches='tight')
plt.close(fig)
print("Saved figures/phase1/fig1_panelD.png/.pdf")
