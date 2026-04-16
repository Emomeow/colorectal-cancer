#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2 Panel B — 沿拟时序分布的关键驱动基因热图
T细胞（上）+ 髓系细胞（下），各含 CRC-fate 基因 + UC-fate 基因
"""

import scanpy as sc
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.sparse import issparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings, io, sys, os
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

plt.rcParams['font.family']  = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

FIG_DIR = 'figures/phase2'
os.makedirs(FIG_DIR, exist_ok=True)


def savefig(fig, name):
    fig.savefig(f'{FIG_DIR}/{name}.png', dpi=200, bbox_inches='tight', facecolor='white')
    fig.savefig(f'{FIG_DIR}/{name}.pdf',           bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  Saved: {FIG_DIR}/{name}.png/.pdf')


# ── 参数 ──────────────────────────────────────────────────────────────
N_BINS   = 50      # 拟时间分 bin 数
SIGMA    = 1.5     # 平滑系数
N_GENES  = 20      # 每个命运取前 N 个基因

# 生物学标记基因（确保在 9009 交集基因集内）
# T细胞 CRC-fate (Tex) marker genes
TEX_MARKERS  = ['HAVCR2', 'ENTPD1', 'TOX', 'TIGIT', 'LAG3', 'PDCD1',
                 'CTLA4', 'CCL5', 'CD7', 'CREM', 'HOPX', 'RGS1']
# T细胞 UC-fate (Tpex/effector) marker genes
TPEX_MARKERS = ['TCF7', 'CCR7', 'SELL', 'LTB', 'JUNB', 'IL7R',
                 'IFITM1', 'S100A4']
# 髓系 CRC-fate (TAM) marker genes
TAM_MARKERS  = ['FTH1', 'SOD2', 'SOCS3', 'OSM', 'HMOX1', 'FCGR1A',
                 'LGALS3', 'CD68', 'HLA-DPB1', 'HLA-DPA1', 'RGS1',
                 'CTSB', 'HLA-DRA', 'HLA-DRB1', 'GPR183']
# 髓系 UC-fate (Mono) marker genes
MONO_MARKERS = ['S100A8', 'S100A9', 'FCN1', 'LYZ', 'CST3',
                 'VCAN', 'CD14', 'TYROBP', 'CTSS']


def get_expr_matrix(adata, genes):
    """提取基因表达矩阵，过滤不在 var_names 中的基因。"""
    valid = [g for g in genes if g in adata.var_names]
    if not valid:
        return None, []
    X = adata[:, valid].X
    if issparse(X):
        X = X.toarray()
    return X, valid


def compute_gene_trends(adata, genes, pseudotime_col='pseudotime',
                         n_bins=N_BINS, sigma=SIGMA):
    """
    按拟时间分 bin，计算各基因平均表达，高斯平滑后返回 (genes × bins) 矩阵。
    """
    pt = adata.obs[pseudotime_col].values
    # 过滤 NaN
    valid_mask = ~np.isnan(pt)
    adata_v = adata[valid_mask]
    pt_v    = pt[valid_mask]

    X, valid_genes = get_expr_matrix(adata_v, genes)
    if X is None:
        return None, []

    # Log-normalize (if raw counts)
    X = np.log1p(X)

    pt_min, pt_max = pt_v.min(), pt_v.max()
    bins = np.linspace(pt_min, pt_max, n_bins + 1)
    bin_idx = np.digitize(pt_v, bins) - 1
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)

    trend = np.zeros((len(valid_genes), n_bins))
    for b in range(n_bins):
        mask = bin_idx == b
        if mask.sum() > 0:
            trend[:, b] = X[mask].mean(axis=0)

    # Gaussian smoothing along pseudotime axis
    for i in range(trend.shape[0]):
        trend[i] = gaussian_filter1d(trend[i], sigma=sigma)

    # Z-score per gene for visualization
    mu  = trend.mean(axis=1, keepdims=True)
    std = trend.std(axis=1, keepdims=True) + 1e-9
    trend_z = (trend - mu) / std

    bin_centers = (bins[:-1] + bins[1:]) / 2
    return trend_z, valid_genes, bin_centers


def find_branch_point(branch_probs, pseudotime, n_bins=N_BINS):
    """
    分支点：fate 概率方差最大处对应的拟时间值。
    """
    pt = pseudotime.values
    valid = ~np.isnan(pt)
    pt_v  = pt[valid]
    bp_v  = branch_probs.values[valid]

    bins    = np.linspace(pt_v.min(), pt_v.max(), n_bins + 1)
    bin_idx = np.clip(np.digitize(pt_v, bins) - 1, 0, n_bins - 1)

    variances = []
    for b in range(n_bins):
        m = bin_idx == b
        if m.sum() > 2:
            variances.append(bp_v[m].var(axis=0).mean())
        else:
            variances.append(0)

    peak_bin = np.argmax(variances)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    return bin_centers[peak_bin]


def plot_heatmap_panel(ax, trend_z, genes, bin_centers, title,
                        branch_pt=None, crc_n=None, uc_n=None,
                        crc_color='#E64B35', uc_color='#4DBBD5'):
    """
    绘制单个拟时间热图。trend_z: (genes × bins)
    """
    vmax = min(np.percentile(np.abs(trend_z), 95), 2.5)
    im = ax.imshow(trend_z, aspect='auto', cmap='RdBu_r',
                   vmin=-vmax, vmax=vmax,
                   extent=[bin_centers[0], bin_centers[-1],
                            len(genes) - 0.5, -0.5])

    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=6.5)
    ax.set_xlabel('Pseudotime', fontsize=8)
    ax.set_title(title, fontsize=9, fontweight='bold', pad=3)
    ax.spines[['top', 'right']].set_visible(False)

    # 分支点竖线
    if branch_pt is not None:
        ax.axvline(branch_pt, color='black', lw=1.2, ls='--', alpha=0.8)
        ax.text(branch_pt, -0.8, 'Branch\npoint',
                ha='center', va='bottom', fontsize=6, color='black')

    # 命运方向标签
    if crc_n is not None:
        ax.text(bin_centers[-1], len(genes) // 4,
                f'CRC-fate\n({crc_n} genes)', ha='right', va='center',
                fontsize=6.5, color=crc_color, fontweight='bold')
    if uc_n is not None:
        ax.text(bin_centers[-1], len(genes) * 3 // 4,
                f'UC-fate\n({uc_n} genes)', ha='right', va='center',
                fontsize=6.5, color=uc_color, fontweight='bold')

    return im


# ══════════════════════════════════════════════════════════════════════
print('Loading data ...')
adata = sc.read_h5ad('processed_data/merged_annotated_fine.h5ad')

# ── T cells ──────────────────────────────────────────────────────────
print('\nComputing T cell gene trends ...')
tcell_mask = adata.obs['cell_type_grouped'] == 'T cells'
adata_t = adata[tcell_mask].copy()

pt_t  = pd.read_csv('phase2_results/tcell_v2/pseudotime_values.csv', index_col=0)
bp_t  = pd.read_csv('phase2_results/tcell_v2/branch_probs.csv', index_col=0)
adata_t.obs['pseudotime'] = pt_t.loc[adata_t.obs_names, '0'].values

# CRC-fate branch 列（Tex方向：SMC07-N对应的列）
# 从数据验证：col 0 = SMC03-T (dominant), col 1 = SMC07-N (alternative)
# 取fate prob variance 最大对应的分支点
branch_pt_t = find_branch_point(bp_t.loc[adata_t.obs_names], adata_t.obs['pseudotime'])
print(f'  T cell branch point pseudotime: {branch_pt_t:.3f}')

# 从 branch DEG 补充额外基因
deg_crc_t = pd.read_csv('phase2_results/tcell_v2/branch_DEG_CRC_fate.csv')
deg_uc_t  = pd.read_csv('phase2_results/tcell_v2/branch_DEG_UC_fate.csv')

# 筛出在 adata.var_names 中的基因
avail = set(adata.var_names)
extra_crc_t = [g for g in deg_crc_t['names'].head(30) if g in avail and g not in TEX_MARKERS][:8]
extra_uc_t  = [g for g in deg_uc_t['names'].head(30)  if g in avail and g not in TPEX_MARKERS][:5]

genes_crc_t = [g for g in TEX_MARKERS  + extra_crc_t if g in avail]
genes_uc_t  = [g for g in TPEX_MARKERS + extra_uc_t  if g in avail]
genes_t     = genes_crc_t + genes_uc_t

trend_t, valid_t, bins_t = compute_gene_trends(adata_t, genes_t)
# Split back into CRC/UC sections
n_crc_t = len([g for g in genes_crc_t if g in valid_t])
n_uc_t  = len([g for g in genes_uc_t  if g in valid_t])
print(f'  T genes: {n_crc_t} CRC-fate + {n_uc_t} UC-fate = {len(valid_t)} total')

# ── Myeloid ──────────────────────────────────────────────────────────
print('\nComputing myeloid gene trends ...')
my_mask = adata.obs['cell_type_grouped'] == 'Myeloids'
adata_m = adata[my_mask].copy()

pt_m  = pd.read_csv('phase2_results/myeloid_v2/pseudotime_values.csv', index_col=0)
bp_m  = pd.read_csv('phase2_results/myeloid_v2/branch_probs.csv', index_col=0)
adata_m.obs['pseudotime'] = pt_m.loc[adata_m.obs_names, '0'].values

branch_pt_m = find_branch_point(bp_m.loc[adata_m.obs_names], adata_m.obs['pseudotime'])
print(f'  Myeloid branch point pseudotime: {branch_pt_m:.3f}')

deg_crc_m = pd.read_csv('phase2_results/myeloid_v2/branch_DEG_CRC_fate.csv')
deg_uc_m  = pd.read_csv('phase2_results/myeloid_v2/branch_DEG_UC_fate.csv')

extra_crc_m = [g for g in deg_crc_m['names'].head(30) if g in avail and g not in TAM_MARKERS][:8]
extra_uc_m  = [g for g in deg_uc_m['names'].head(30)  if g in avail and g not in MONO_MARKERS][:5]

genes_crc_m = [g for g in TAM_MARKERS  + extra_crc_m if g in avail]
genes_uc_m  = [g for g in MONO_MARKERS + extra_uc_m  if g in avail]
genes_m     = genes_crc_m + genes_uc_m

trend_m, valid_m, bins_m = compute_gene_trends(adata_m, genes_m)
n_crc_m = len([g for g in genes_crc_m if g in valid_m])
n_uc_m  = len([g for g in genes_uc_m  if g in valid_m])
print(f'  Myeloid genes: {n_crc_m} CRC-fate + {n_uc_m} UC-fate = {len(valid_m)} total')


# ══════════════════════════════════════════════════════════════════════
#  Figure
# ══════════════════════════════════════════════════════════════════════
print('\nPlotting Panel B ...')

fig = plt.figure(figsize=(14, 10))
gs  = fig.add_gridspec(2, 2, hspace=0.4, wspace=0.35,
                        left=0.12, right=0.97, top=0.93, bottom=0.06,
                        width_ratios=[1, 0.04])

# Row 0: T cells
ax_t  = fig.add_subplot(gs[0, 0])
ax_tc = fig.add_subplot(gs[0, 1])   # colorbar

im_t = plot_heatmap_panel(
    ax_t, trend_t, valid_t, bins_t,
    title='T cell Gene Dynamics along Pseudotime',
    branch_pt=branch_pt_t,
    crc_n=n_crc_t, uc_n=n_uc_t,
)
plt.colorbar(im_t, cax=ax_tc, label='Z-score')
ax_tc.tick_params(labelsize=6)
ax_tc.yaxis.label.set_fontsize(7)

# 在 CRC/UC 基因分界处画横线
if trend_t is not None and n_crc_t > 0 and n_uc_t > 0:
    ax_t.axhline(n_crc_t - 0.5, color='black', lw=0.8, ls='-', alpha=0.5)
    ax_t.text(bins_t[0], n_crc_t // 2, '← CRC-fate genes',
              va='center', ha='left', fontsize=6.5, color='#E64B35',
              fontweight='bold')
    ax_t.text(bins_t[0], n_crc_t + n_uc_t // 2, '← UC-fate genes',
              va='center', ha='left', fontsize=6.5, color='#4DBBD5',
              fontweight='bold')

# Row 1: Myeloid
ax_m  = fig.add_subplot(gs[1, 0])
ax_mc = fig.add_subplot(gs[1, 1])

im_m = plot_heatmap_panel(
    ax_m, trend_m, valid_m, bins_m,
    title='Myeloid Gene Dynamics along Pseudotime',
    branch_pt=branch_pt_m,
    crc_n=n_crc_m, uc_n=n_uc_m,
)
plt.colorbar(im_m, cax=ax_mc, label='Z-score')
ax_mc.tick_params(labelsize=6)
ax_mc.yaxis.label.set_fontsize(7)

if trend_m is not None and n_crc_m > 0 and n_uc_m > 0:
    ax_m.axhline(n_crc_m - 0.5, color='black', lw=0.8, ls='-', alpha=0.5)
    ax_m.text(bins_m[0], n_crc_m // 2, '← CRC-fate genes',
              va='center', ha='left', fontsize=6.5, color='#E64B35',
              fontweight='bold')
    ax_m.text(bins_m[0], n_crc_m + n_uc_m // 2, '← UC-fate genes',
              va='center', ha='left', fontsize=6.5, color='#4DBBD5',
              fontweight='bold')

fig.text(0.005, 0.97, 'B', fontsize=14, fontweight='bold', va='top')

savefig(fig, 'figure2_panelB')
print('Done.')
