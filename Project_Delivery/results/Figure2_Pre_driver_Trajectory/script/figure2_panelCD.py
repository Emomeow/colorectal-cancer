#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2 Panel C + D — GSE144735 空间定位模拟
Panel C: UMAP 展示 Border/Tumor/Normal 分区 + Pre-driver 评分
Panel D: Infiltration Score 条形图（各区 Pre-driver 亚群丰度）
"""

import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import warnings, io, sys, os
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

plt.rcParams['font.family']  = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

FIG_DIR = 'figures/phase2'
RES_DIR = 'phase2_results/spatial'
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)


def savefig(fig, name):
    fig.savefig(f'{FIG_DIR}/{name}.png', dpi=200, bbox_inches='tight', facecolor='white')
    fig.savefig(f'{FIG_DIR}/{name}.pdf',           bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  Saved: {FIG_DIR}/{name}.png/.pdf')


# ── 配色 ──────────────────────────────────────────────────────────────
CLASS_COLORS = {
    'Normal': '#4DBBD5',
    'Border': '#F39B7F',
    'Tumor':  '#E64B35',
}
CELLTYPE_COLORS = {
    'T cells':          '#3C5488',
    'Myeloids':         '#E64B35',
    'B cells':          '#00A087',
    'Epithelial cells': '#7E6148',
    'Stromal cells':    '#8491B4',
    'Mast cells':       '#B09C85',
}
CLASS_ORDER = ['Normal', 'Border', 'Tumor']

# ══════════════════════════════════════════════════════════════════════
print('Loading GSE144735 ...')
adata = sc.read_h5ad('processed_data/GSE144735.h5ad')
print(f'  {adata.n_obs:,} cells, Class: {adata.obs["Class"].value_counts().to_dict()}')

# ── UMAP ──────────────────────────────────────────────────────────────
print('Computing UMAP ...')
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat_v3')
sc.pp.pca(adata, n_comps=30, use_highly_variable=True)
sc.pp.neighbors(adata, n_neighbors=20, use_rep='X_pca')
sc.tl.umap(adata, min_dist=0.3, random_state=42)
xy = adata.obsm['X_umap']

# ══════════════════════════════════════════════════════════════════════
#  Panel C — UMAP (2×2 subplot)
# ══════════════════════════════════════════════════════════════════════
print('Plotting Panel C ...')
fig_c, axes = plt.subplots(2, 2, figsize=(12, 10))
fig_c.subplots_adjust(hspace=0.3, wspace=0.25,
                       left=0.06, right=0.97, top=0.93, bottom=0.06)

# C1 — 按 Class 着色
ax = axes[0, 0]
for cls in CLASS_ORDER:
    m = adata.obs['Class'] == cls
    ax.scatter(xy[m, 0], xy[m, 1], c=CLASS_COLORS[cls],
               s=1.5, alpha=0.6, linewidths=0, rasterized=True, label=cls)
patches = [mpatches.Patch(color=CLASS_COLORS[c], label=c) for c in CLASS_ORDER]
ax.legend(handles=patches, fontsize=7, frameon=False, markerscale=2, loc='upper right')
ax.set_title('Spatial Region', fontsize=10, fontweight='bold')
ax.set_xticks([]); ax.set_yticks([])
ax.spines[['top','right','bottom','left']].set_visible(False)

# C2 — 按 Cell_type 着色
ax = axes[0, 1]
for ct in list(CELLTYPE_COLORS.keys()):
    m = adata.obs['Cell_type'] == ct
    ax.scatter(xy[m, 0], xy[m, 1], c=CELLTYPE_COLORS[ct],
               s=1.5, alpha=0.6, linewidths=0, rasterized=True, label=ct)
patches2 = [mpatches.Patch(color=CELLTYPE_COLORS[c], label=c) for c in CELLTYPE_COLORS]
ax.legend(handles=patches2, fontsize=6.5, frameon=False, loc='upper right')
ax.set_title('Cell Type', fontsize=10, fontweight='bold')
ax.set_xticks([]); ax.set_yticks([])
ax.spines[['top','right','bottom','left']].set_visible(False)

# C3 — TAM_score 热度图
ax = axes[1, 0]
tam_vals = adata.obs['TAM_score'].values
vmax = np.percentile(tam_vals, 98)
vmin = np.percentile(tam_vals, 2)
sc3 = ax.scatter(xy[:, 0], xy[:, 1], c=tam_vals, cmap='YlOrRd',
                  vmin=vmin, vmax=vmax, s=1.5, alpha=0.7,
                  linewidths=0, rasterized=True)
plt.colorbar(sc3, ax=ax, shrink=0.7, pad=0.02, label='TAM score')
ax.set_title('TAM Signature Score', fontsize=10, fontweight='bold')
ax.set_xticks([]); ax.set_yticks([])
ax.spines[['top','right','bottom','left']].set_visible(False)

# C4 — CD8Tex_score 热度图
ax = axes[1, 1]
tex_vals = adata.obs['CD8Tex_score'].values
vmax2 = np.percentile(tex_vals, 98)
vmin2 = np.percentile(tex_vals, 2)
sc4 = ax.scatter(xy[:, 0], xy[:, 1], c=tex_vals, cmap='PuRd',
                  vmin=vmin2, vmax=vmax2, s=1.5, alpha=0.7,
                  linewidths=0, rasterized=True)
plt.colorbar(sc4, ax=ax, shrink=0.7, pad=0.02, label='CD8_Tex score')
ax.set_title('CD8_Tex Exhaustion Score', fontsize=10, fontweight='bold')
ax.set_xticks([]); ax.set_yticks([])
ax.spines[['top','right','bottom','left']].set_visible(False)

fig_c.text(0.005, 0.97, 'C', fontsize=14, fontweight='bold', va='top')
savefig(fig_c, 'figure2_panelC')


# ══════════════════════════════════════════════════════════════════════
#  Panel D — Infiltration Score 条形图
#  定义: 各 patient × class 中 Myeloid / T cell 的 Pre-driver 亚群比例
# ══════════════════════════════════════════════════════════════════════
print('Computing Infiltration Scores ...')

def infiltration_score(adata_sub, score_col, threshold_pct=75):
    """
    每个 patient 内各区的 Pre-driver 比例：
    细胞评分 > 第 threshold_pct 百分位数 视为 Pre-driver
    """
    threshold = np.percentile(adata_sub.obs[score_col], threshold_pct)
    adata_sub = adata_sub.copy()
    adata_sub.obs['is_predriver'] = adata_sub.obs[score_col] > threshold

    records = []
    for patient in adata_sub.obs['Patient'].unique():
        for cls in CLASS_ORDER:
            m = (adata_sub.obs['Patient'] == patient) & (adata_sub.obs['Class'] == cls)
            if m.sum() < 10:
                continue
            frac = adata_sub.obs.loc[m, 'is_predriver'].mean()
            records.append({'Patient': patient, 'Class': cls, 'fraction': frac})
    return pd.DataFrame(records)


# ── TAM Infiltration（仅 Myeloid 细胞）────────────────────────────────
adata_my = adata[adata.obs['Cell_type'] == 'Myeloids'].copy()
df_tam = infiltration_score(adata_my, 'TAM_score', threshold_pct=75)
df_tam['module'] = 'TAM'

# ── CD8_Tex Infiltration（仅 T 细胞）────────────────────────────────
adata_t = adata[adata.obs['Cell_type'] == 'T cells'].copy()
df_tex = infiltration_score(adata_t, 'CD8Tex_score', threshold_pct=75)
df_tex['module'] = 'CD8_Tex'

df_all = pd.concat([df_tam, df_tex], ignore_index=True)
df_all.to_csv(f'{RES_DIR}/infiltration_scores.csv', index=False)

# 统计检验: Border vs Tumor
print('\nMann-Whitney U tests:')
results_stats = []
for module in ['TAM', 'CD8_Tex']:
    sub = df_all[df_all['module'] == module]
    for c1, c2 in [('Normal', 'Border'), ('Border', 'Tumor'), ('Normal', 'Tumor')]:
        v1 = sub[sub['Class'] == c1]['fraction'].values
        v2 = sub[sub['Class'] == c2]['fraction'].values
        if len(v1) >= 3 and len(v2) >= 3:
            stat, p = mannwhitneyu(v1, v2, alternative='two-sided')
            results_stats.append({'module': module, 'comparison': f'{c1} vs {c2}', 'p': p})
            print(f'  {module}  {c1} vs {c2}: p={p:.4f}')

# ── 作图 ─────────────────────────────────────────────────────────────
print('\nPlotting Panel D ...')

fig_d, axes_d = plt.subplots(1, 2, figsize=(10, 5))
fig_d.subplots_adjust(wspace=0.35, left=0.1, right=0.97, top=0.90, bottom=0.12)

module_titles = {'TAM': 'TAM Infiltration Score\n(Myeloid cells)',
                 'CD8_Tex': 'CD8_Tex Infiltration Score\n(T cells)'}
module_colors = {'TAM': '#E64B35', 'CD8_Tex': '#3C5488'}

for idx, module in enumerate(['TAM', 'CD8_Tex']):
    ax = axes_d[idx]
    sub = df_all[df_all['module'] == module]

    # 各区均值 + 误差
    means  = sub.groupby('Class')['fraction'].mean().reindex(CLASS_ORDER)
    sems   = sub.groupby('Class')['fraction'].sem().reindex(CLASS_ORDER)
    n_pts  = sub.groupby('Class')['fraction'].count().reindex(CLASS_ORDER)

    x = np.arange(len(CLASS_ORDER))
    bars = ax.bar(x, means, color=[CLASS_COLORS[c] for c in CLASS_ORDER],
                  alpha=0.85, width=0.55, edgecolor='white', linewidth=0.5,
                  yerr=sems, capsize=4, error_kw={'linewidth': 1.2})

    # 散点（per patient）
    for ci, cls in enumerate(CLASS_ORDER):
        pts = sub[sub['Class'] == cls]['fraction'].values
        jitter = np.random.default_rng(42).uniform(-0.12, 0.12, len(pts))
        ax.scatter(ci + jitter, pts, c='black', s=14, alpha=0.6,
                   linewidths=0, zorder=3)

    # 显著性标注
    stat_sub = [r for r in results_stats if r['module'] == module]
    y_max = means.max() + sems.max() + 0.02
    sig_pairs = [('Normal', 'Border'), ('Border', 'Tumor'), ('Normal', 'Tumor')]
    y_offsets = [0.02, 0.06, 0.10]
    for (c1, c2), y_off in zip(sig_pairs, y_offsets):
        r = next((r for r in stat_sub if r['comparison'] == f'{c1} vs {c2}'), None)
        if r is None: continue
        p = r['p']
        if p < 0.001:  sig = '***'
        elif p < 0.01: sig = '**'
        elif p < 0.05: sig = '*'
        else:          sig = 'ns'

        x1 = CLASS_ORDER.index(c1)
        x2 = CLASS_ORDER.index(c2)
        y  = y_max + y_off
        ax.plot([x1, x1, x2, x2], [y-0.005, y, y, y-0.005],
                color='black', lw=0.8)
        ax.text((x1+x2)/2, y+0.003, sig, ha='center', va='bottom', fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(CLASS_ORDER, fontsize=9)
    ax.set_ylabel('Pre-driver Fraction\n(cells above 75th percentile)', fontsize=8)
    ax.set_title(module_titles[module], fontsize=9, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)

    # 患者数标注
    for ci, cls in enumerate(CLASS_ORDER):
        ax.text(ci, -0.018, f'n={int(n_pts[cls])}',
                ha='center', va='top', fontsize=6.5, color='grey')

fig_d.text(0.005, 0.97, 'D', fontsize=14, fontweight='bold', va='top')
savefig(fig_d, 'figure2_panelD')

print('\nAll done.')
