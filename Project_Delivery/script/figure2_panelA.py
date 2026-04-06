#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2 Panel A — Palantir 拟时序演化轨迹图
T细胞 + 髓系细胞  各3个子图：亚型 | 拟时间 | 命运分支概率
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import warnings, io, sys, os
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

plt.rcParams['font.family']  = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

FIG_DIR = 'figures/phase2'
os.makedirs(FIG_DIR, exist_ok=True)

# ── 配色 ──────────────────────────────────────────────────────────────
TCELL_COLORS = {
    'CD8_Tex':    '#D62728',  # tab:red
    'CD8_Tpex':   '#FF7F0E',  # tab:orange
    'CD8_Teff':   '#1F77B4',  # tab:blue
    'Treg':       '#2CA02C',  # tab:green
    'Th17':       '#9467BD',  # tab:purple
    'CD4_Tconv':  '#8C564B',  # tab:brown
    'Tfh':        '#E377C2',  # tab:pink
    'gdT':        '#17BECF',  # tab:cyan
    'Unknown_T':  '#CCCCCC',
}
MYELOID_COLORS = {
    'TAM':            '#D62728',  # tab:red
    'Inflam_Mac':     '#FF7F0E',  # tab:orange
    'Classical_Mono': '#1F77B4',  # tab:blue
    'NonClass_Mono':  '#2CA02C',  # tab:green
    'cDC1':           '#9467BD',  # tab:purple
    'cDC2':           '#17BECF',  # tab:cyan
    'pDC':            '#E377C2',  # tab:pink
    'Unknown_My':     '#CCCCCC',
}

PSEUDOTIME_CMAP = 'viridis'
FATE_CMAP       = 'RdBu_r'


def savefig(fig, name):
    fig.savefig(f'{FIG_DIR}/{name}.png', dpi=200, bbox_inches='tight', facecolor='white')
    fig.savefig(f'{FIG_DIR}/{name}.pdf',           bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  Saved: {FIG_DIR}/{name}.png/.pdf')


def compute_umap(adata_sub, n_neighbors=30, min_dist=0.3, seed=42):
    """用 X_scVI 给子集计算专属 UMAP，保留更细的内部结构。"""
    sc.pp.neighbors(adata_sub, use_rep='X_scvi', n_neighbors=n_neighbors)
    sc.tl.umap(adata_sub, min_dist=min_dist, random_state=seed)
    return adata_sub


def scatter_umap(ax, xy, color_vals, cmap, vmin=None, vmax=None,
                 s=2, alpha=0.5, rasterized=True):
    """通用 UMAP 散点（连续色）。"""
    sc_kw = dict(s=s, alpha=alpha, linewidths=0, rasterized=rasterized)
    if vmin is None: vmin = np.nanpercentile(color_vals, 2)
    if vmax is None: vmax = np.nanpercentile(color_vals, 98)
    sc = ax.scatter(xy[:, 0], xy[:, 1], c=color_vals,
                    cmap=cmap, vmin=vmin, vmax=vmax, **sc_kw)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[['top','right','bottom','left']].set_visible(False)
    return sc


def scatter_umap_cat(ax, xy, labels, color_map, s=2, alpha=0.6):
    """通用 UMAP 散点（分类色），顺序：先画灰色未知，再画有色亚型。"""
    cats = list(color_map.keys())
    # grey first
    for cat in [c for c in cats if 'Unknown' in c or 'Mast' in c]:
        m = labels == cat
        if m.any():
            ax.scatter(xy[m, 0], xy[m, 1], c=color_map[cat],
                       s=s, alpha=0.3, linewidths=0, rasterized=True)
    # colored
    for cat in [c for c in cats if 'Unknown' not in c and 'Mast' not in c]:
        m = labels == cat
        if m.any():
            ax.scatter(xy[m, 0], xy[m, 1], c=color_map[cat],
                       s=s, alpha=alpha, linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[['top','right','bottom','left']].set_visible(False)


def add_legend(ax, color_map, present_cats, fontsize=6, ncol=1):
    handles = [Line2D([0], [0], marker='o', color='w',
                      markerfacecolor=color_map[c], markersize=5, label=c)
               for c in color_map if c in present_cats]
    ax.legend(handles=handles, fontsize=fontsize, frameon=False,
              loc='upper left', ncol=ncol,
              handlelength=0.8, handleheight=0.8, borderpad=0.3)


def add_colorbar(fig, ax, sc_obj, label, fontsize=7):
    cbar = fig.colorbar(sc_obj, ax=ax, shrink=0.7, pad=0.02, aspect=20)
    cbar.set_label(label, fontsize=fontsize)
    cbar.ax.tick_params(labelsize=6)


# ══════════════════════════════════════════════════════════════════════
print('Loading data ...')
adata = sc.read_h5ad('processed_data/merged_annotated_fine.h5ad')
print(f'  {adata.n_obs:,} cells')

# ── T cell subset ────────────────────────────────────────────────────
print('\nPreparing T cell trajectory ...')
tcell_mask = adata.obs['cell_type_grouped'] == 'T cells'
adata_t = adata[tcell_mask].copy()

# Join Palantir results
pt_t  = pd.read_csv('phase2_results/tcell_v2/pseudotime_values.csv', index_col=0)
bp_t  = pd.read_csv('phase2_results/tcell_v2/branch_probs.csv', index_col=0)

adata_t.obs['pseudotime'] = pt_t.loc[adata_t.obs_names, '0']

# Identify CRC-fate column (terminal state in CRC tissue with highest CRC bias)
# From analysis: col 0 (SMC03-T, GSE132465 tumor) = dominant fate
# col 1 (SMC07-N, GSE132465 normal) = alternative fate
# Label CRC_fate = col that corresponds to CD8_Tex enrichment
# (SMC03-T is tumor; in myeloid the CRC fate was SMC171-T-SING)
# Both T cell terminals are GSE132465(CRC); use the one with higher prob in CD8_Tex cells
tex_mask = adata_t.obs['cell_subtype_fine'] == 'CD8_Tex'
col0_tex = bp_t.loc[adata_t.obs_names[tex_mask], bp_t.columns[0]].mean()
col1_tex = bp_t.loc[adata_t.obs_names[tex_mask], bp_t.columns[1]].mean()
crc_col_t = bp_t.columns[0] if col0_tex >= col1_tex else bp_t.columns[1]
print(f'  T cell CRC-fate column (Tex-enriched): {crc_col_t[:40]}')

adata_t.obs['fate_CRC'] = bp_t.loc[adata_t.obs_names, crc_col_t].values

# UMAP
print('  Computing T cell sub-UMAP ...')
adata_t = compute_umap(adata_t)
xy_t = adata_t.obsm['X_umap']
subtypes_t = adata_t.obs['cell_subtype_fine'].values

# ── Myeloid subset ────────────────────────────────────────────────────
print('\nPreparing myeloid trajectory ...')
my_mask = adata.obs['cell_type_grouped'] == 'Myeloids'
adata_m = adata[my_mask].copy()

pt_m  = pd.read_csv('phase2_results/myeloid_v2/pseudotime_values.csv', index_col=0)
bp_m  = pd.read_csv('phase2_results/myeloid_v2/branch_probs.csv', index_col=0)

adata_m.obs['pseudotime'] = pt_m.loc[adata_m.obs_names, '0']

# CRC-fate column: SMC171-T-SING (GSE132257, CRC tumor) — highest CRC bias
crc_col_m = [c for c in bp_m.columns if 'GSE132257' in c or 'GSE132465' in c][0]
print(f'  Myeloid CRC-fate column: {crc_col_m[:45]}')
adata_m.obs['fate_CRC'] = bp_m.loc[adata_m.obs_names, crc_col_m].values

# UMAP
print('  Computing myeloid sub-UMAP ...')
adata_m = compute_umap(adata_m, n_neighbors=20)
xy_m = adata_m.obsm['X_umap']
subtypes_m = adata_m.obs['cell_subtype_fine'].values


# ══════════════════════════════════════════════════════════════════════
#  Figure: 2 rows × 3 cols
# ══════════════════════════════════════════════════════════════════════
print('\nPlotting Panel A ...')

fig = plt.figure(figsize=(14, 9))
gs  = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.12,
                        left=0.04, right=0.97, top=0.93, bottom=0.04)

row_labels = ['T cells', 'Myeloids']
col_titles = ['Cell Subtype', 'Pseudotime', 'Suppressive Fate Probability']

# ── Row 0: T cells ─────────────────────────────────────────────────
ax00 = fig.add_subplot(gs[0, 0])
ax01 = fig.add_subplot(gs[0, 1])
ax02 = fig.add_subplot(gs[0, 2])

# Col 0 — subtype
present_t = np.unique(subtypes_t)
cmap_t_use = {k: v for k, v in TCELL_COLORS.items() if k in present_t}
scatter_umap_cat(ax00, xy_t, subtypes_t, cmap_t_use, s=1.5, alpha=0.7)
add_legend(ax00, cmap_t_use, present_t, fontsize=6)

# Col 1 — pseudotime
sc01 = scatter_umap(ax01, xy_t, adata_t.obs['pseudotime'].values,
                    cmap=PSEUDOTIME_CMAP, s=1.5)
add_colorbar(fig, ax01, sc01, 'Pseudotime')

# Col 2 — fate probability
sc02 = scatter_umap(ax02, xy_t, adata_t.obs['fate_CRC'].values,
                    cmap='YlOrRd', vmin=0, vmax=1, s=1.5)
add_colorbar(fig, ax02, sc02, 'CD8_Tex fate prob.')

# ── Row 1: Myeloids ────────────────────────────────────────────────
ax10 = fig.add_subplot(gs[1, 0])
ax11 = fig.add_subplot(gs[1, 1])
ax12 = fig.add_subplot(gs[1, 2])

present_m = np.unique(subtypes_m)
cmap_m_use = {k: v for k, v in MYELOID_COLORS.items() if k in present_m}
scatter_umap_cat(ax10, xy_m, subtypes_m, cmap_m_use, s=3, alpha=0.8)
add_legend(ax10, cmap_m_use, present_m, fontsize=6)

sc11 = scatter_umap(ax11, xy_m, adata_m.obs['pseudotime'].values,
                    cmap=PSEUDOTIME_CMAP, s=3)
add_colorbar(fig, ax11, sc11, 'Pseudotime')

sc12 = scatter_umap(ax12, xy_m, adata_m.obs['fate_CRC'].values,
                    cmap='YlOrRd', vmin=0, vmax=1, s=3)
add_colorbar(fig, ax12, sc12, 'TAM fate prob.')

# ── Labels ──────────────────────────────────────────────────────────
for col_i, title in enumerate(col_titles):
    fig.text(0.05 + col_i * 0.315, 0.955, title,
             ha='center', va='top', fontsize=10, fontweight='bold')

for row_i, label in enumerate(row_labels):
    y = 0.91 - row_i * 0.46
    fig.text(0.003, y, label, ha='left', va='center',
             fontsize=9, fontweight='bold', rotation=90)

# Panel label
fig.text(0.005, 0.97, 'A', fontsize=14, fontweight='bold', va='top')

savefig(fig, 'figure2_panelA')
print('\nDone.')
