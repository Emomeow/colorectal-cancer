#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2 Panel E — TF Regulon 活性分析（decoupleR + CollecTRI）
在 TAM / CD8_Tex 亚群中推断 TF 活性，识别 Pre-driver 核心转录调控网络
"""

import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import issparse
import decoupler as dc
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
RES_DIR = 'phase2_results/regulon'
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)


def savefig(fig, name):
    fig.savefig(f'{FIG_DIR}/{name}.png', dpi=200, bbox_inches='tight', facecolor='white')
    fig.savefig(f'{FIG_DIR}/{name}.pdf',           bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  Saved: {FIG_DIR}/{name}.png/.pdf')


# ── 关注的 TF（免疫耗竭 / 肿瘤 TAM 相关）────────────────────────────
TF_FOCUS = [
    # 耗竭/抑制
    'TOX', 'NR4A1', 'NFATC1', 'BATF', 'IRF4',
    # 炎症 & 免疫抑制
    'STAT3', 'HIF1A', 'RELA', 'NFKB1', 'CEBPB',
    # TAM 极化
    'IRF1', 'STAT1', 'PPARG', 'KLF4', 'MYC',
    # AP-1 家族
    'JUN', 'FOS', 'FOSL1', 'FOSL2',
    # 额外
    'SP1', 'TP53', 'RUNX1',
]

# 目标亚群及对比亚群
TARGET_SUBTYPES = ['TAM', 'CD8_Tex']
COMPARE_TO = {
    'TAM':     'Classical_Mono',   # TAM vs 单核细胞
    'CD8_Tex': 'CD8_Teff',         # 耗竭 vs 效应 T
}

# ══════════════════════════════════════════════════════════════════════
print('Loading data ...')
adata = sc.read_h5ad('processed_data/merged_annotated_fine.h5ad')
print(f'  {adata.n_obs:,} cells')

# ── 获取 CollecTRI TF-target 网络 ────────────────────────────────────
print('Loading CollecTRI network ...')
net = dc.get_collectri(organism='human', split_complexes=False)
# 仅保留 target 在 adata.var_names 中的记录
avail_genes = set(adata.var_names)
net = net[net['target'].isin(avail_genes)]
# 仅保留 focus TF（以及有足够靶基因的 TF）
tf_counts = net.groupby('source')['target'].count()
valid_tfs  = tf_counts[tf_counts >= 5].index
net = net[net['source'].isin(valid_tfs)]
print(f'  Network: {len(net):,} TF-target pairs, {net["source"].nunique()} TFs')

# ── 子集：仅取目标亚群 + 对比亚群 ─────────────────────────────────────
subtypes_needed = set(TARGET_SUBTYPES + list(COMPARE_TO.values()))
mask = adata.obs['cell_subtype_fine'].isin(subtypes_needed)
adata_sub = adata[mask].copy()
print(f'  Subset: {adata_sub.n_obs:,} cells ({adata_sub.obs["cell_subtype_fine"].value_counts().to_dict()})')

# 准备 log-normalized 表达矩阵
sc.pp.normalize_total(adata_sub, target_sum=1e4)
sc.pp.log1p(adata_sub)

# ── 运行 ULM（单变量线性模型）推断 TF 活性 ───────────────────────────
print('Running ULM TF activity inference ...')
dc.run_ulm(
    mat   = adata_sub,
    net   = net,
    source= 'source',
    target= 'target',
    weight= 'weight',
    verbose=False,
    use_raw=False,
)
# ULM 结果存在 obsm['ulm_estimate'] 和 obsm['ulm_pvals']
acts = adata_sub.obsm['ulm_estimate']   # DataFrame: cells × TFs
print(f'  TF activity matrix: {acts.shape}')

# ── 筛出 focus TF（取交集）────────────────────────────────────────────
focus_avail = [tf for tf in TF_FOCUS if tf in acts.columns]
print(f'  Focus TFs available: {len(focus_avail)}/{len(TF_FOCUS)}  → {focus_avail}')
acts_focus = acts[focus_avail]

# ── 保存 TF 活性均值 per subtype ────────────────────────────────────
acts_focus_df = acts_focus.copy()
acts_focus_df['subtype'] = adata_sub.obs['cell_subtype_fine'].values
mean_by_subtype = acts_focus_df.groupby('subtype')[focus_avail].mean()
mean_by_subtype.to_csv(f'{RES_DIR}/tf_activity_mean.csv')
print('\nMean TF activity (selected):')
print(mean_by_subtype.T.round(2).to_string())

# ══════════════════════════════════════════════════════════════════════
#  Plot Panel E — 2 部分
#  E1: 热图 — TAM vs Mono 的 TF 活性差值
#  E2: 热图 — CD8_Tex vs CD8_Teff 的 TF 活性差值
# ══════════════════════════════════════════════════════════════════════
print('\nPlotting Panel E ...')

# 每对比计算 delta activity (target - compare)
def compute_delta(acts_df, subtype_col, target, compare, tfs):
    t_mean = acts_df[acts_df[subtype_col] == target][tfs].mean()
    c_mean = acts_df[acts_df[subtype_col] == compare][tfs].mean()
    delta  = t_mean - c_mean
    return delta.sort_values(ascending=False)


acts_focus_df2 = acts_focus.copy()
acts_focus_df2['subtype'] = adata_sub.obs['cell_subtype_fine'].values

delta_tam = compute_delta(acts_focus_df2, 'subtype', 'TAM',     'Classical_Mono', focus_avail)
delta_tex = compute_delta(acts_focus_df2, 'subtype', 'CD8_Tex', 'CD8_Teff',        focus_avail)
delta_tam.to_csv(f'{RES_DIR}/delta_tf_TAM_vs_Mono.csv', header=['delta'])
delta_tex.to_csv(f'{RES_DIR}/delta_tf_Tex_vs_Teff.csv', header=['delta'])

# ─── 组合热图：focus TF × {TAM-Mono, Tex-Teff} ───────────────────────
# 按 TAM delta 排序，显示一致性
tfs_plot = delta_tam.index.tolist()   # 按TAM delta排序
mat_plot  = pd.DataFrame({
    'TAM vs\nClassical_Mono': delta_tam[tfs_plot],
    'CD8_Tex vs\nCD8_Teff':   delta_tex[tfs_plot],
})

fig, axes = plt.subplots(1, 2, figsize=(11, 8),
                          gridspec_kw={'width_ratios': [2, 0.05],
                                        'wspace': 0.04})
ax    = axes[0]
ax_cb = axes[1]

vmax = np.percentile(np.abs(mat_plot.values), 95)
im   = ax.imshow(mat_plot.values, aspect='auto', cmap='RdBu_r',
                  vmin=-vmax, vmax=vmax)

ax.set_xticks(range(len(mat_plot.columns)))
ax.set_xticklabels(mat_plot.columns, fontsize=9)
ax.set_yticks(range(len(tfs_plot)))
ax.set_yticklabels(tfs_plot, fontsize=8)
ax.set_title('TF Regulon Activity\n(Pre-driver vs Control subtype)',
              fontsize=10, fontweight='bold', pad=6)

# 标注正负区域背景
for i, tf in enumerate(tfs_plot):
    v_tam = mat_plot.loc[tf, 'TAM vs\nClassical_Mono']
    v_tex = mat_plot.loc[tf, 'CD8_Tex vs\nCD8_Teff']
    # 在两列右侧标注最大值方向
    max_v = max(abs(v_tam), abs(v_tex))
    if max_v > vmax * 0.5:
        color = '#B22222' if (v_tam + v_tex) > 0 else '#1A6B8A'
        ax.text(1.55, i, '◀', fontsize=7, va='center', color=color, fontweight='bold')

plt.colorbar(im, cax=ax_cb, label='ΔActivity\n(z-score)')
ax_cb.tick_params(labelsize=7)
ax_cb.yaxis.label.set_fontsize(8)

# 分区域说明
ax.axvline(0.5, color='white', lw=2)  # 分隔线

# 高亮关键 TF 名称
key_tfs = {'STAT3', 'HIF1A', 'CEBPB', 'IRF1', 'TOX', 'NR4A1', 'BATF'}
for i, tf in enumerate(tfs_plot):
    if tf in key_tfs:
        ax.get_yticklabels()[i].set_fontweight('bold')
        ax.get_yticklabels()[i].set_color('#B22222')

fig.text(0.005, 0.97, 'E', fontsize=14, fontweight='bold', va='top')

# 添加小型注释说明
fig.text(0.08, 0.02,
         'Red = higher activity in Pre-driver subtype  |  '
         'Blue = lower activity in Pre-driver subtype  |  '
         'Network: CollecTRI (ULM)',
         fontsize=7, ha='left', color='grey', style='italic')

fig.subplots_adjust(left=0.18, right=0.93, top=0.93, bottom=0.08)
savefig(fig, 'figure2_panelE')
print('Done.')
