#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Root Robustness Check — Phase 2 Pseudotime Sensitivity Analysis
===============================================================
验证拟时序轨迹方向对根细胞选择的稳健性。

当前使用的 root 策略：UC 细胞中 DC1 最小值（扩散空间最"原始"端）。
本脚本额外测试两种中性标准：
  1. 最低细胞周期评分（S + G2M 分）→ 最静止 / 最未分化
  2. 最少总 UMI（最简单基因表达）

对每个 root 独立运行 Palantir，输出：
  - 拟时序排名 Spearman 相关矩阵
  - UMAP 可视化（3 个 root 方案并排）
  - branch point 伪时序值比较表

输出：Project_Delivery/results/Figure2_Pre_driver_Trajectory/Trajectory_Analysis/figures/
      root_robustness_*.png / .pdf
"""

import scanpy as sc
import numpy as np
import pandas as pd
import palantir
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import warnings, io, sys, os
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

plt.rcParams['font.family']  = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

FIG_DIR = ('Project_Delivery/results/Figure2_Pre_driver_Trajectory/'
           'Trajectory_Analysis/figures')
os.makedirs(FIG_DIR, exist_ok=True)

CELL_TYPE = 'Myeloids'   # 以髓系为示例（T 细胞可将此改为 'T cells'）
BATCHES   = ['GSE125527', 'GSE132465', 'GSE132257', 'GSE150115']

# Cell cycle genes (subset; used for S/G2M scoring)
S_GENES   = ['MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG',
              'GINS2','MCM6','CDCA7','DTL','PRIM1','UHRF1','MLF1IP',
              'HELLS','RFC2','RPA2','NASP','RAD51AP1']
G2M_GENES = ['HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A',
              'NDC80','CKS2','NUF2','CKS1B','MKI67','TMPO','CENPF',
              'TACC3','FAM64A','SMC4','CCNB2','CKAP2L','CKAP2']


def savefig(fig, name):
    fig.savefig(f'{FIG_DIR}/{name}.png', dpi=200, bbox_inches='tight', facecolor='white')
    fig.savefig(f'{FIG_DIR}/{name}.pdf',           bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  Saved: {FIG_DIR}/{name}.png/.pdf')


def run_palantir_from_root(ms_data, root_cell, n_waypoints=800):
    """Run Palantir and return pseudotime Series."""
    pr_res = palantir.core.run_palantir(
        ms_data, early_cell=root_cell,
        num_waypoints=n_waypoints,
        use_early_cell_as_start=True,
        n_jobs=4,
    )
    return pr_res.pseudotime, pr_res.branch_probs


def find_branch_point(branch_probs, pseudotime, n_bins=50):
    """Pseudotime at maximum fate-probability variance (branch point)."""
    pt = pseudotime.values
    valid = ~np.isnan(pt)
    pt_v, bp_v = pt[valid], branch_probs.values[valid]
    bins    = np.linspace(pt_v.min(), pt_v.max(), n_bins + 1)
    bin_idx = np.clip(np.digitize(pt_v, bins) - 1, 0, n_bins - 1)
    variances = [bp_v[bin_idx == b].var(axis=0).mean()
                 if (bin_idx == b).sum() > 2 else 0
                 for b in range(n_bins)]
    peak = np.argmax(variances)
    return ((bins[:-1] + bins[1:]) / 2)[peak]


# ── Load data ─────────────────────────────────────────────────────────
print('Loading data ...')
adata_all = sc.read_h5ad('processed_data/merged_annotated_fine.h5ad')

mask = (adata_all.obs['cell_type_grouped'] == CELL_TYPE) & \
       (adata_all.obs['batch'].isin(BATCHES))
adata = adata_all[mask].copy()
print(f'  {CELL_TYPE}: {adata.n_obs:,} cells')

# ── Normalize for cell cycle scoring ─────────────────────────────────
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.score_genes_cell_cycle(adata, s_genes=S_GENES, g2m_genes=G2M_GENES)
adata.obs['cc_score'] = adata.obs['S_score'] + adata.obs['G2M_score']
print(f'  Cell cycle score: {adata.obs["cc_score"].describe().to_dict()}')

# ── Build diffusion map (shared across all root tests) ────────────────
print('Computing diffusion map ...')
dm_res = palantir.utils.run_diffusion_maps(
    pd.DataFrame(adata.obsm['X_scvi'], index=adata.obs_names),
    n_components=10, knn=30,
)
ms_data = palantir.utils.determine_multiscale_space(dm_res)
dc1 = dm_res['EigenVectors'].iloc[:, 0]

# UMAP for visualization
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=30)
sc.tl.umap(adata, min_dist=0.3, random_state=0)
umap = adata.obsm['X_umap']

# ── Define 3 root candidates ──────────────────────────────────────────
uc_mask = adata.obs['disease'] == 'UC'

# Root 1 (original): UC cell at min DC1
root_dc1 = dc1[uc_mask].idxmin()
print(f'\nRoot 1 (DC1 min):       {root_dc1[:50]}')

# Root 2: UC cell with lowest cell cycle score
root_cc = adata.obs.loc[uc_mask, 'cc_score'].idxmin()
print(f'Root 2 (min cc_score):  {root_cc[:50]}')

# Root 3: UC cell with lowest total UMI (n_counts or total_counts)
umi_col = 'total_counts' if 'total_counts' in adata.obs.columns else 'n_counts'
if umi_col in adata.obs.columns:
    root_umi = adata.obs.loc[uc_mask, umi_col].idxmin()
else:
    # fallback: use sum of raw counts if available
    raw_sum = np.array(adata.X.sum(axis=1)).ravel()
    adata.obs['_raw_sum'] = raw_sum
    root_umi = adata.obs.loc[uc_mask, '_raw_sum'].idxmin()
print(f'Root 3 (min UMI):       {root_umi[:50]}')

ROOTS = {
    'Root1_DC1min':       root_dc1,
    'Root2_minCellCycle': root_cc,
    'Root3_minUMI':       root_umi,
}

# ── Run Palantir for each root ────────────────────────────────────────
results = {}
for name, root in ROOTS.items():
    print(f'\nRunning Palantir — {name} ...')
    pt, bp = run_palantir_from_root(ms_data, root)
    results[name] = {'pseudotime': pt, 'branch_probs': bp}
    bp_val = find_branch_point(bp, pt)
    print(f'  Branch point pseudotime: {bp_val:.3f}')
    print(f'  Pseudotime range: {pt.min():.3f} – {pt.max():.3f}')

# ── Spearman rank correlation between pseudotimes ─────────────────────
print('\n── Pseudotime rank correlations ──')
pt_df = pd.DataFrame({n: r['pseudotime'] for n, r in results.items()})
corr_matrix = pt_df.corr(method='spearman')
print(corr_matrix.to_string())

# ── Plot: 3×3 UMAP grid (root × pseudotime/batch/branch_pt) ──────────
fig, axes = plt.subplots(3, 3, figsize=(14, 13))
fig.suptitle(f'Root Robustness Check — {CELL_TYPE}\n'
             f'Spearman rank correlations between pseudotimes',
             fontsize=11, fontweight='bold')

shuf = np.random.default_rng(0).permutation(adata.n_obs)

for row_i, (name, res) in enumerate(results.items()):
    pt_vals = res['pseudotime'].reindex(adata.obs_names).values
    root_id = list(ROOTS.values())[row_i]
    root_idx = adata.obs_names.get_loc(root_id)

    # Col 0: pseudotime
    ax = axes[row_i, 0]
    sc_obj = ax.scatter(umap[shuf, 0], umap[shuf, 1],
                        c=pt_vals[shuf], s=1.5, alpha=0.5,
                        cmap='viridis', linewidths=0)
    ax.scatter(umap[root_idx, 0], umap[root_idx, 1],
               s=200, c='red', marker='*', zorder=10)
    plt.colorbar(sc_obj, ax=ax, shrink=0.7, pad=0.02)
    ax.set_title(f'{name}\nPseudotime', fontsize=8, fontweight='bold')
    ax.axis('off')

    # Col 1: subtype
    ax = axes[row_i, 1]
    COLORS = {'TAM': '#D62728', 'Inflam_Mac': '#FF7F0E',
              'Classical_Mono': '#1F77B4', 'NonClass_Mono': '#2CA02C',
              'cDC1': '#9467BD', 'cDC2': '#17BECF', 'pDC': '#E377C2',
              'CD8_Tex': '#D62728', 'CD8_Tpex': '#FF7F0E',
              'CD8_Teff': '#1F77B4', 'Treg': '#2CA02C',
              'Th17': '#9467BD', 'CD4_Tconv': '#8C564B',
              'Tfh': '#E377C2', 'gdT': '#17BECF'}
    subtypes = adata.obs.get('cell_subtype_fine',
                              adata.obs.get('cell_type_grouped')).values
    for sub, col in COLORS.items():
        m = subtypes == sub
        if m.any():
            ax.scatter(umap[m, 0], umap[m, 1], c=col,
                       s=1.5, alpha=0.6, linewidths=0, label=sub)
    ax.set_title('Cell subtype', fontsize=8, fontweight='bold')
    ax.axis('off')
    if row_i == 0:
        ax.legend(fontsize=5, markerscale=3, frameon=False,
                  loc='upper left', ncol=2)

    # Col 2: branch probability (first terminal)
    ax = axes[row_i, 2]
    bp_vals = res['branch_probs'].iloc[:, 0].reindex(adata.obs_names).values
    bp_pt   = find_branch_point(res['branch_probs'], res['pseudotime'])
    sc_bp = ax.scatter(umap[shuf, 0], umap[shuf, 1],
                       c=bp_vals[shuf], s=1.5, alpha=0.5,
                       cmap='YlOrRd', vmin=0, vmax=1, linewidths=0)
    plt.colorbar(sc_bp, ax=ax, shrink=0.7, pad=0.02)
    ax.set_title(f'Branch prob → Term1\nBP={bp_pt:.3f}', fontsize=8, fontweight='bold')
    ax.axis('off')

# Spearman table as text
corr_str = (f"Spearman ρ:\n"
            f"Root1 vs Root2: {corr_matrix.iloc[0,1]:.3f}\n"
            f"Root1 vs Root3: {corr_matrix.iloc[0,2]:.3f}\n"
            f"Root2 vs Root3: {corr_matrix.iloc[1,2]:.3f}")
fig.text(0.5, 0.01, corr_str, ha='center', va='bottom',
         fontsize=9, fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout(rect=[0, 0.06, 1, 1])
savefig(fig, 'root_robustness_umap')

# ── Summary table ─────────────────────────────────────────────────────
summary_rows = []
for name, res in results.items():
    pt  = res['pseudotime']
    bp  = find_branch_point(res['branch_probs'], pt)
    summary_rows.append({
        'Root strategy': name,
        'Root cell (truncated)': list(ROOTS.values())[list(results.keys()).index(name)][:40],
        'Pseudotime min': f'{pt.min():.3f}',
        'Pseudotime max': f'{pt.max():.3f}',
        'Branch point pt': f'{bp:.3f}',
        'Spearman vs Root1': f'{corr_matrix.loc[name, "Root1_DC1min"]:.3f}',
    })
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(f'{FIG_DIR}/root_robustness_summary.csv', index=False)
print('\n── Summary ──')
print(summary_df.to_string(index=False))

print(f'\nConclusion:')
rho12 = corr_matrix.iloc[0, 1]
rho13 = corr_matrix.iloc[0, 2]
if min(rho12, rho13) > 0.8:
    print(f'  Trajectory direction ROBUST (min ρ={min(rho12,rho13):.3f} > 0.8)')
elif min(rho12, rho13) > 0.6:
    print(f'  Trajectory largely consistent (min ρ={min(rho12,rho13):.3f}); minor sensitivity to root')
else:
    print(f'  Trajectory SENSITIVE to root selection (min ρ={min(rho12,rho13):.3f}); consider alternative strategy')

print('\nDone.')
