#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplementary Figure 2 — 特异性亚群 Marker 基因扩充列表
• 用 Wilcoxon rank_genes_groups 计算数据驱动 marker
• 与文献 curated 基因融合，取 top 10 per subtype
• Dotplot: 行=subtype, 列=基因（按亚群分组）
• 点大小 = fraction expressing, 颜色 = mean expression z-score
"""

import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import issparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import warnings, io, sys, os
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

plt.rcParams['font.family']  = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

FIG_DIR = 'figures/phase2'
RES_DIR = 'phase2_results/markers'
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

N_MARKERS = 10   # markers per subtype in dotplot

# ── 文献 curated 核心 marker（优先级高）────────────────────────────
CURATED = {
    'TAM':            ['APOE', 'SPP1', 'FTH1', 'CTSB', 'SOD2', 'HMOX1',
                       'LGALS3', 'CD68', 'HLA-DRA', 'FCGR1A', 'RGS1', 'SOCS3'],
    'Inflam_Mac':     ['IL1B', 'CXCL8', 'IL6', 'CCL2', 'TNF',
                       'CXCL3', 'G0S2', 'PTGS2', 'NFKBIA', 'VEGFA'],
    'Classical_Mono': ['S100A8', 'S100A9', 'FCN1', 'LYZ', 'CST3',
                       'VCAN', 'CD14', 'TYROBP', 'CTSS', 'MNDA'],
    'cDC2':           ['FCER1A', 'CD1C', 'CLEC10A', 'CLNK',
                       'PKIB', 'CD86', 'HLA-DQA1', 'SPI1'],
    'pDC':            ['IRF7', 'TCF4', 'PLD4', 'PLAC8',
                       'IRF8', 'GZMB', 'TSPAN13', 'ITM2C'],
    'CD8_Tex':        ['HAVCR2', 'ENTPD1', 'PDCD1', 'TIGIT', 'LAG3',
                       'CTLA4', 'TOX', 'CREM', 'HOPX', 'CXCL13'],
    'CD8_Tpex':       ['TCF7', 'CCR7', 'SELL', 'IL7R', 'LTB',
                       'TOX2', 'CXCR5', 'BCL6', 'SLAMF6', 'JUNB'],
    'CD8_Teff':       ['GZMB', 'PRF1', 'IFNG', 'NKG7', 'CCL4',
                       'FGFBP2', 'CX3CR1', 'TBX21', 'KLRG1', 'GZMA'],
    'CD4_Tconv':      ['IL7R', 'LTB', 'MAL', 'LDHB', 'CCR7',
                       'NOSIP', 'TCF7', 'AQP3'],
    'Treg':           ['FOXP3', 'IL2RA', 'CTLA4', 'TNFRSF9', 'TIGIT',
                       'IKZF2', 'RTKN2', 'BATF', 'LAYN', 'DUSP4'],
    'Th17':           ['IL17A', 'RORC', 'CCR6', 'RORA',
                       'CAPG', 'FURIN', 'KLRB1'],
    'Tfh':            ['CXCR5', 'BCL6', 'ICOS', 'PDCD1', 'TOX2',
                       'CD200', 'BTLA', 'MAF', 'SH2D1A'],
    'gdT':            ['TRDC', 'TRGC2', 'GNLY', 'KLRB1',
                       'NKG7', 'FCGR3A', 'GZMB', 'CXCR6'],
    'NK cells':       ['GNLY', 'NKG7', 'KLRB1', 'NCAM1', 'FCGR3A',
                       'TYROBP', 'FGFBP2', 'PRF1', 'KLRD1'],
    'B cells':        ['MS4A1', 'CD79A', 'CD79B', 'IGHM', 'IGHD',
                       'BANK1', 'RALGPS2', 'FCRL1', 'CD22', 'PAX5'],
    'Mast cells':     ['SRGN', 'CD63', 'ANXA1', 'FCER1G', 'ALOX5AP',
                       'LAPTM4A', 'CD9', 'SAMSN1', 'TYROBP'],
    'Stromal':        ['IGFBP7', 'LGALS1', 'COL6A2', 'GSN', 'TIMP1',
                       'IGFBP4', 'VIM', 'IFITM3', 'HSPB1'],
    'Epithelial':     ['KRT18', 'S100A6', 'ATP1B1', 'TPM1', 'SPINT2',
                       'KLF5', 'MGST1', 'PERP', 'CYSTM1'],
    'Malignant':      ['KRT18', 'SOX4', 'GNAS', 'LGALS3', 'TXN',
                       'PERP', 'KLF5', 'LLGL2', 'PLCG2', 'SLC12A2'],
}

LINEAGE_GROUPS = {
    'Myeloid':        ['TAM', 'Inflam_Mac', 'Classical_Mono', 'cDC2', 'pDC'],
    'T cells':        ['CD8_Tex', 'CD8_Tpex', 'CD8_Teff', 'CD4_Tconv',
                       'Treg', 'Th17', 'Tfh', 'gdT'],
    'NK / B / Other': ['NK cells', 'B cells', 'Mast cells',
                       'Stromal', 'Epithelial', 'Malignant'],
}
LINEAGE_COLORS = {
    'Myeloid':        '#E64B35',
    'T cells':        '#3C5488',
    'NK / B / Other': '#00A087',
}

# ══════════════════════════════════════════════════════════════════════
print('Loading data ...')
adata = sc.read_h5ad('processed_data/merged_annotated_fine.h5ad')
print(f'  {adata.n_obs:,} cells, {adata.n_vars:,} genes')

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

avail_genes = set(adata.var_names)
avail_subtypes = [s for g in LINEAGE_GROUPS.values() for s in g
                  if s in adata.obs['cell_subtype_fine'].values]

# ── rank_genes_groups（Wilcoxon）────────────────────────────────────
print('Running Wilcoxon rank_genes_groups ...')
adata_sub = adata[adata.obs['cell_subtype_fine'].isin(avail_subtypes)].copy()
sc.tl.rank_genes_groups(adata_sub, groupby='cell_subtype_fine',
                         method='wilcoxon', n_genes=30, key_added='markers')
print('  Done.')

# ── 融合 curated + computed，构建 final marker list per subtype ─────
print('Building final marker sets ...')
final_markers = {}
marker_table_rows = []

for st in avail_subtypes:
    df = sc.get.rank_genes_groups_df(adata_sub, group=st, key='markers')
    # computed top markers (significant, logFC > 0)
    computed = df[(df['pvals_adj'] < 0.05) & (df['logfoldchanges'] > 0)]['names'].tolist()

    # curated (only those in var_names)
    curated_avail = [g for g in CURATED.get(st, []) if g in avail_genes]

    # merge: curated first, then fill with computed (no duplicates)
    merged = list(curated_avail)
    for g in computed:
        if g not in merged:
            merged.append(g)
    final_markers[st] = merged[:N_MARKERS]

    n_cur = len(curated_avail)
    n_comp = len([g for g in final_markers[st] if g not in curated_avail])
    print(f'  {st}: {n_cur} curated + {n_comp} computed = {len(final_markers[st])} total')

    for g in final_markers[st]:
        row = df[df['names'] == g]
        lfc = float(row['logfoldchanges'].values[0]) if len(row) else np.nan
        padj = float(row['pvals_adj'].values[0]) if len(row) else np.nan
        marker_table_rows.append({'subtype': st, 'gene': g,
                                   'logFC': round(lfc, 3), 'padj': padj,
                                   'source': 'curated' if g in curated_avail else 'computed'})

# Save marker table
marker_df = pd.DataFrame(marker_table_rows)
marker_df.to_csv(f'{RES_DIR}/subtype_markers.csv', index=False)
print(f'  Marker table saved → {RES_DIR}/subtype_markers.csv')

# ── 预计算 mean expression + fraction expressing ────────────────────
print('Computing expression statistics ...')
all_genes_flat = list(dict.fromkeys(
    g for st in avail_subtypes for g in final_markers[st]
))

mean_mat = pd.DataFrame(index=avail_subtypes, columns=all_genes_flat, dtype=float)
frac_mat = pd.DataFrame(index=avail_subtypes, columns=all_genes_flat, dtype=float)

for st in avail_subtypes:
    mask = adata.obs['cell_subtype_fine'] == st
    X = adata[mask, all_genes_flat].X
    if issparse(X): X = X.toarray()
    mean_mat.loc[st] = X.mean(axis=0)
    frac_mat.loc[st] = (X > 0).mean(axis=0)

# Z-score per gene (across all subtypes)
z_mat = mean_mat.copy()
for g in all_genes_flat:
    col = mean_mat[g]
    z_mat[g] = (col - col.mean()) / (col.std() + 1e-9)
z_mat = z_mat.clip(-2.5, 2.5)

# ══════════════════════════════════════════════════════════════════════
#  Dotplot function
# ══════════════════════════════════════════════════════════════════════
def dotplot_lineage(ax, lineage_name, subtypes, lineage_color,
                     max_dot_size=280, cmap='RdBu_r', vmax=2.5):
    gene_list   = []
    gene_bounds = {}   # st -> (start, end) col indices
    col = 0
    for st in subtypes:
        g = [x for x in final_markers[st] if x in all_genes_flat]
        gene_bounds[st] = (col, col + len(g) - 1)
        for gg in g:
            if gg not in gene_list:
                gene_list.append(gg)
        col += len(g)

    n_genes = len(gene_list)
    n_sub   = len(subtypes)

    xs, ys, ss, cs = [], [], [], []
    for si, st in enumerate(subtypes):
        for gi, gene in enumerate(gene_list):
            frac = float(frac_mat.loc[st, gene]) if gene in frac_mat.columns else 0
            zval = float(z_mat.loc[st, gene])    if gene in z_mat.columns    else 0
            xs.append(gi); ys.append(si)
            ss.append(frac * max_dot_size)
            cs.append(zval)

    sc_obj = ax.scatter(xs, ys, s=ss, c=cs,
                        cmap=cmap, vmin=-vmax, vmax=vmax,
                        linewidths=0.2, edgecolors='#bbbbbb', zorder=3)

    ax.set_xlim(-0.6, n_genes - 0.4)
    ax.set_ylim(-0.6, n_sub  - 0.4)

    ax.set_xticks(range(n_genes))
    ax.set_xticklabels(gene_list, rotation=55, ha='right', fontsize=7.8)
    ax.set_yticks(range(n_sub))
    n_cells = {st: (adata.obs['cell_subtype_fine'] == st).sum() for st in subtypes}
    ax.set_yticklabels(
        [f'{st}  (n={n_cells[st]:,})' for st in subtypes],
        fontsize=9.0
    )
    ax.tick_params(axis='x', length=2, pad=1)
    ax.tick_params(axis='y', length=0)

    # Alternating band per subtype gene block + subtype label at top
    for ki, st in enumerate(subtypes):
        g0, g1 = gene_bounds[st]
        alpha = 0.10 if ki % 2 == 0 else 0.03
        ax.axvspan(g0 - 0.5, g1 + 0.5, alpha=alpha,
                   color=lineage_color, zorder=0)
        mid = (g0 + g1) / 2
        ax.text(mid, n_sub - 0.35, st,
                ha='center', va='bottom', fontsize=7.8,
                color=lineage_color, fontweight='bold', rotation=30,
                clip_on=False)
        if ki > 0:
            ax.axvline(g0 - 0.5, color='#cccccc', lw=0.7, zorder=1)

    for si in range(n_sub - 1):
        ax.axhline(si + 0.5, color='#eeeeee', lw=0.5, zorder=1)

    ax.set_title(f'  {lineage_name}', fontsize=11, fontweight='bold',
                  color=lineage_color, loc='left', pad=30)
    ax.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
    ax.set_facecolor('white')

    return sc_obj


# ── 图形尺寸计算 ────────────────────────────────────────────────────
n_my  = len(LINEAGE_GROUPS['Myeloid'])
n_t   = len(LINEAGE_GROUPS['T cells'])
n_oth = len(LINEAGE_GROUPS['NK / B / Other'])

ROW_H  = 0.46   # inches per subtype row
HEADER = 1.6    # header space per panel (for rotated gene labels)
GAP    = 0.5    # gap between panels
YPAD   = 0.8    # bottom/top padding

fig_h = (n_my + n_t + n_oth) * ROW_H + 3 * HEADER + 2 * GAP + YPAD
fig_w = 25.0
CB_W  = 0.18   # colorbar width inches

print('Plotting ...')
fig = plt.figure(figsize=(fig_w, fig_h))
fig.patch.set_facecolor('white')

# ── layout via axes with manual position ────────────────────────────
left   = 0.14
right  = 0.955
cb_l   = 0.960
cb_w   = CB_W / fig_w

def panel_bottom_height(n_rows):
    return n_rows * ROW_H / fig_h, (n_rows * ROW_H + HEADER) / fig_h

h_oth_plot, h_oth_tot = panel_bottom_height(n_oth)
h_t_plot,   h_t_tot   = panel_bottom_height(n_t)
h_my_plot,  h_my_tot  = panel_bottom_height(n_my)

bot_oth = YPAD / fig_h / 2
bot_t   = bot_oth + h_oth_tot + GAP / fig_h
bot_my  = bot_t   + h_t_tot   + GAP / fig_h

ax_my  = fig.add_axes([left, bot_my  + (h_my_tot  - h_my_plot),  right - left, h_my_plot])
ax_t   = fig.add_axes([left, bot_t   + (h_t_tot   - h_t_plot),   right - left, h_t_plot])
ax_oth = fig.add_axes([left, bot_oth + (h_oth_tot  - h_oth_plot), right - left, h_oth_plot])
ax_cb  = fig.add_axes([cb_l,  bot_oth, cb_w,
                        bot_my + h_my_tot - bot_oth])

sc1 = dotplot_lineage(ax_my,  'Myeloid',
                       LINEAGE_GROUPS['Myeloid'],        LINEAGE_COLORS['Myeloid'])
sc2 = dotplot_lineage(ax_t,   'T cells',
                       LINEAGE_GROUPS['T cells'],        LINEAGE_COLORS['T cells'])
sc3 = dotplot_lineage(ax_oth, 'NK / B / Other',
                       LINEAGE_GROUPS['NK / B / Other'], LINEAGE_COLORS['NK / B / Other'])

# Colorbar
cb = plt.colorbar(sc1, cax=ax_cb)
cb.set_label('Mean Expression\n(z-score, clipped ±2.5)', fontsize=8.5, labelpad=6)
ax_cb.tick_params(labelsize=8)

# Dot size legend
fracs_demo  = [0.1, 0.25, 0.5, 0.75, 1.0]
labels_demo = ['10%', '25%', '50%', '75%', '100%']
legend_handles = [
    Line2D([0], [0], marker='o', color='w',
           markerfacecolor='#444444', markeredgecolor='#aaaaaa',
           markersize=np.sqrt(f * 280) * 0.75,
           label=lbl)
    for f, lbl in zip(fracs_demo, labels_demo)
]
ax_my.legend(
    handles=legend_handles,
    title='Fraction\nexpressing',
    title_fontsize=8.5,
    fontsize=8,
    loc='upper left',
    bbox_to_anchor=(1.01, 1.0),
    frameon=True,
    edgecolor='#cccccc',
    borderpad=0.7,
    handlelength=0.6,
)

# Global title + subtitle
fig.text(0.5, 0.998,
         'Supplementary Figure 2 — Cell Subtype-Specific Marker Genes',
         fontsize=13.5, fontweight='bold', ha='center', va='top')
fig.text(0.5, 0.993,
         'Dot size: fraction of cells expressing gene  |  '
         'Color: mean log-normalized expression (z-score per gene, ±2.5 clip)  |  '
         'Markers: literature-curated + Wilcoxon top DEGs',
         fontsize=8.5, ha='center', va='top', color='grey', style='italic')

out_name = 'supp_fig2_markers'
fig.savefig(f'{FIG_DIR}/{out_name}.png', dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(f'{FIG_DIR}/{out_name}.pdf',           bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {FIG_DIR}/{out_name}.png/.pdf')
print('Done.')
