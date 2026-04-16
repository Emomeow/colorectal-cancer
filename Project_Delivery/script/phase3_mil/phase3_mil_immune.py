#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 3 — 免疫细胞专属 ABMIL
Bags 仅含 TNKILC + Myeloid 细胞（排除上皮/基质），
强制模型从免疫微环境特征分类 MMRp vs MMRd。
核心问题：在 MMRp bags 中，哪个免疫亚群的注意力权重最高？
"""

import scanpy as sc
import numpy as np
import pandas as pd
import io, sys, os, warnings
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import Adam
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

torch.manual_seed(42)
np.random.seed(42)

plt.rcParams['font.family']  = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

FIG_DIR = 'figures/phase3'
RES_DIR = 'phase3_results/mil'
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

# ── 超参数 ─────────────────────────────────────────────────────────
N_PCA      = 50
MAX_CELLS  = 1500
HIDDEN_DIM = 256
ATT_DIM    = 128
DROPOUT    = 0.25
LR         = 1e-4
WD         = 1e-5
EPOCHS     = 100
N_FOLDS    = 5
PATIENCE   = 20

# ══════════════════════════════════════════════════════════════════════
#  Step 1: 加载 joint 数据，筛选免疫细胞
# ══════════════════════════════════════════════════════════════════════
print('Loading joint data ...')
adata_full = sc.read_h5ad('processed_data/GSE178341_joint.h5ad')
print(f'  {adata_full.n_obs:,} cells')

IMMUNE_TYPES = ['TNKILC', 'Myeloid']
adata = adata_full[adata_full.obs['clTopLevel'].isin(IMMUNE_TYPES)].copy()
print(f'  Immune cells: {adata.n_obs:,}')
print(f'  MMRp: {(adata.obs["MMRStatus"]=="MMRp").sum():,}  '
      f'MMRd: {(adata.obs["MMRStatus"]=="MMRd").sum():,}')

# 过滤患者：要求每患者至少 50 个免疫细胞
pid_counts = adata.obs.groupby(['PID', 'MMRStatus']).size()
valid_pids = pid_counts[pid_counts >= 50].index.get_level_values('PID').unique()
adata = adata[adata.obs['PID'].isin(valid_pids)].copy()
print(f'  After per-patient filter (≥50 cells): {adata.n_obs:,} cells, '
      f'{adata.obs["PID"].nunique()} patients')

# 重新计算 HVG + PCA（仅免疫细胞）
print('  Recomputing HVG + PCA on immune cells only ...')
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3',
                              layer='counts', subset=False)
sc.pp.pca(adata, n_comps=N_PCA, use_highly_variable=True)
# UMAP
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=40)
sc.tl.umap(adata, min_dist=0.3, random_state=42)
print('  Done.')

# ── Bag 构建 ────────────────────────────────────────────────────────
pca_feat = adata.obsm['X_pca']
pids     = adata.obs['PID'].values
mmr_lbl  = (adata.obs['MMRStatus'] == 'MMRd').astype(int).values

unique_pids = adata.obs['PID'].unique()
patient_bags = {}
for pid in unique_pids:
    idx = np.where(pids == pid)[0]
    lbl = mmr_lbl[idx[0]]
    patient_bags[pid] = {'idx': idx, 'label': lbl}

all_pids   = np.array(list(patient_bags.keys()))
all_labels = np.array([patient_bags[p]['label'] for p in all_pids])
n_mmrd = all_labels.sum()
n_mmrp = (1 - all_labels).sum()
print(f'\nBags: {len(all_pids)} (MMRd={n_mmrd}, MMRp={n_mmrp})')
print(f'Cells/bag: median={np.median([len(patient_bags[p]["idx"]) for p in all_pids]):.0f}')

# ══════════════════════════════════════════════════════════════════════
#  ABMIL Model
# ══════════════════════════════════════════════════════════════════════
class GatedABMIL(nn.Module):
    def __init__(self, input_dim=N_PCA, hidden_dim=HIDDEN_DIM,
                 att_dim=ATT_DIM, dropout=DROPOUT):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim), nn.ReLU(),
        )
        self.att_V = nn.Sequential(nn.Linear(hidden_dim, att_dim), nn.Tanh())
        self.att_U = nn.Sequential(nn.Linear(hidden_dim, att_dim), nn.Sigmoid())
        self.att_w = nn.Linear(att_dim, 1, bias=False)
        self.classifier = nn.Sequential(nn.Dropout(dropout),
                                         nn.Linear(hidden_dim, 1))

    def forward(self, x, return_attention=False):
        H = self.encoder(x)
        A = self.att_w(self.att_V(H) * self.att_U(H))
        A = F.softmax(A, dim=0)
        M = (A * H).sum(dim=0, keepdim=True)
        logit = self.classifier(M).squeeze()
        if return_attention:
            return logit, A.squeeze()
        return logit


def train_epoch(model, optimizer, pids_list, pca_feat, patient_bags,
                max_cells=MAX_CELLS):
    model.train()
    total_loss = 0
    order = pids_list.copy()
    np.random.shuffle(order)
    for pid in order:
        idx = patient_bags[pid]['idx']
        label = torch.tensor(patient_bags[pid]['label'], dtype=torch.float32)
        if len(idx) > max_cells:
            idx = np.random.choice(idx, max_cells, replace=False)
        x = torch.tensor(pca_feat[idx], dtype=torch.float32)
        optimizer.zero_grad()
        loss = F.binary_cross_entropy_with_logits(model(x), label)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
    return total_loss / len(order)


@torch.no_grad()
def eval_bags(model, pids_eval, pca_feat, patient_bags):
    model.eval()
    probs, labels = [], []
    for pid in pids_eval:
        idx = patient_bags[pid]['idx']
        x   = torch.tensor(pca_feat[idx], dtype=torch.float32)
        probs.append(torch.sigmoid(model(x)).item())
        labels.append(patient_bags[pid]['label'])
    return np.array(probs), np.array(labels)


# ══════════════════════════════════════════════════════════════════════
#  5-Fold CV
# ══════════════════════════════════════════════════════════════════════
print('\n5-Fold CV (immune-only bags) ...')
skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)
fold_aucs   = []
best_states = {}

for fold, (tr_idx, va_idx) in enumerate(skf.split(all_pids, all_labels)):
    pids_tr = all_pids[tr_idx]
    pids_va = all_pids[va_idx]
    model = GatedABMIL()
    opt   = Adam(model.parameters(), lr=LR, weight_decay=WD)
    best_auc, best_state, patience_cnt = 0, None, 0
    for epoch in range(EPOCHS):
        train_epoch(model, opt, list(pids_tr), pca_feat, patient_bags)
        probs_va, lbl_va = eval_bags(model, pids_va, pca_feat, patient_bags)
        auc = roc_auc_score(lbl_va, probs_va) \
              if len(np.unique(lbl_va)) > 1 else 0.5
        if auc > best_auc:
            best_auc = auc
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            patience_cnt = 0
        else:
            patience_cnt += 1
        if patience_cnt >= PATIENCE:
            break
    print(f'  Fold {fold+1}: best AUC={best_auc:.4f} (epoch {epoch+1})')
    fold_aucs.append(best_auc)
    best_states[fold] = best_state

print(f'  Mean AUC: {np.mean(fold_aucs):.4f} ± {np.std(fold_aucs):.4f}')
pd.DataFrame({'fold': range(1,6), 'auc': fold_aucs}).to_csv(
    f'{RES_DIR}/cv_auc_immune.csv', index=False)

# ══════════════════════════════════════════════════════════════════════
#  最终模型 + 全量注意力提取
# ══════════════════════════════════════════════════════════════════════
print('\nFinal model + attention extraction ...')
best_fold = int(np.argmax(fold_aucs))
final_model = GatedABMIL()
final_model.load_state_dict(best_states[best_fold])
# Fine-tune 30 epochs on full data
opt_f = Adam(final_model.parameters(), lr=LR*0.5, weight_decay=WD)
for _ in range(30):
    train_epoch(final_model, opt_f, list(all_pids), pca_feat, patient_bags)
torch.save(final_model.state_dict(), f'{RES_DIR}/abmil_immune_final.pt')

# Extract attention
final_model.eval()
att_scores = np.zeros(adata.n_obs, dtype=np.float32)
with torch.no_grad():
    for pid in all_pids:
        idx = patient_bags[pid]['idx']
        x   = torch.tensor(pca_feat[idx], dtype=torch.float32)
        _, A = final_model(x, return_attention=True)
        att_scores[idx] = A.numpy()

adata.obs['attention_score'] = att_scores

# ══════════════════════════════════════════════════════════════════════
#  统计分析
# ══════════════════════════════════════════════════════════════════════
print('\nStatistical analysis ...')

# 按亚群统计注意力分值（分 MMRp/MMRd）
def subtype_att_stats(adata_sub, mmr):
    sub = adata_sub.obs[adata_sub.obs['MMRStatus'] == mmr]
    g = sub.groupby('clMidwayPr')['attention_score']
    df = g.agg(['mean','median','std','count']).reset_index()
    df.columns = ['subtype','mean_att','median_att','std_att','n_cells']
    return df[df['n_cells'] >= 50].sort_values('mean_att', ascending=False)

stats_mmrp = subtype_att_stats(adata, 'MMRp')
stats_mmrd = subtype_att_stats(adata, 'MMRd')
print('\nMMRp top attention subtypes:')
print(stats_mmrp.head(10).to_string(index=False))
print('\nMMRd top attention subtypes:')
print(stats_mmrd.head(10).to_string(index=False))

stats_mmrp.to_csv(f'{RES_DIR}/immune_att_MMRp.csv', index=False)
stats_mmrd.to_csv(f'{RES_DIR}/immune_att_MMRd.csv', index=False)

# Pre-driver 亚群 vs all others (MMRp 患者)
mmrp_obs = adata.obs[adata.obs['MMRStatus'] == 'MMRp']
macro_att  = mmrp_obs[mmrp_obs['clMidwayPr'] == 'Macro']['attention_score'].values
tcd8_att   = mmrp_obs[mmrp_obs['clMidwayPr'] == 'TCD8']['attention_score'].values
other_att  = mmrp_obs[~mmrp_obs['clMidwayPr'].isin(['Macro','TCD8'])]['attention_score'].values

tests = {}
for name, group_att in [('Macro', macro_att), ('TCD8', tcd8_att)]:
    if len(group_att) > 0:
        stat, p = mannwhitneyu(group_att, other_att, alternative='greater')
        tests[name] = p
        print(f'  MMRp — {name} vs others: p={p:.4e}')

# 亚群间两两比较 (top 免疫亚群, MMRp)
subtypes_test = ['Macro', 'TCD8', 'TCD4', 'Mono', 'DC', 'NK', 'Tgd']
pairs_results = []
for st in subtypes_test:
    grp = mmrp_obs[mmrp_obs['clMidwayPr'] == st]['attention_score'].values
    rest = mmrp_obs[mmrp_obs['clMidwayPr'] != st]['attention_score'].values
    if len(grp) >= 30:
        _, p = mannwhitneyu(grp, rest, alternative='two-sided')
        pairs_results.append({'subtype': st, 'n': len(grp), 'p_vs_rest': p,
                               'mean_att': grp.mean()})
pairs_df = pd.DataFrame(pairs_results)
if len(pairs_df) > 0:
    _, pairs_df['fdr'], _, _ = multipletests(pairs_df['p_vs_rest'], method='fdr_bh')
    pairs_df = pairs_df.sort_values('mean_att', ascending=False)
    print('\nAll pairwise tests (MMRp bags):')
    print(pairs_df.to_string(index=False))
    pairs_df.to_csv(f'{RES_DIR}/immune_att_pairwise_MMRp.csv', index=False)

# ══════════════════════════════════════════════════════════════════════
#  Figure 3B — 免疫专属 MIL 注意力分析
# ══════════════════════════════════════════════════════════════════════
print('\nPlotting Figure 3B ...')

SUBTYPE_COLORS = {
    'Macro':   '#E64B35', 'TCD8':  '#3C5488', 'Mono':    '#F39B7F',
    'DC':      '#E64B35', 'TCD4':  '#3C5488', 'NK':      '#00A087',
    'B':       '#00A087', 'Tgd':   '#8491B4', 'TZBTB16': '#8491B4',
    'Mast':    '#B09C85', 'ILC':   '#7E6148', 'Granulo': '#B09C85',
}
MMR_COLORS = {'MMRp': '#E64B35', 'MMRd': '#4DBBD5'}

xy = adata.obsm['X_umap']

fig = plt.figure(figsize=(20, 16))
fig.patch.set_facecolor('white')
gs = gridspec.GridSpec(3, 3, figure=fig,
                       hspace=0.42, wspace=0.35,
                       left=0.07, right=0.97, top=0.92, bottom=0.06)

# 1. UMAP — cell subtype
ax = fig.add_subplot(gs[0, 0])
subtypes_uniq = adata.obs['clMidwayPr'].unique()
for st in subtypes_uniq:
    m = adata.obs['clMidwayPr'] == st
    col = SUBTYPE_COLORS.get(st, '#aaaaaa')
    ax.scatter(xy[m, 0], xy[m, 1], c=col, s=1.2, alpha=0.5,
               linewidths=0, rasterized=True, label=st)
ax.legend(fontsize=6.5, frameon=False, markerscale=3, loc='upper right',
           ncol=2, handletextpad=0.3)
ax.set_title('Immune Cell Subtypes (clMidwayPr)', fontsize=10, fontweight='bold')
ax.set_xticks([]); ax.set_yticks([])
ax.spines[['top','right','bottom','left']].set_visible(False)

# 2. UMAP — MMR status
ax = fig.add_subplot(gs[0, 1])
for mmr, color in MMR_COLORS.items():
    m = adata.obs['MMRStatus'] == mmr
    ax.scatter(xy[m, 0], xy[m, 1], c=color, s=1.2, alpha=0.5,
               linewidths=0, rasterized=True, label=mmr)
ax.legend(fontsize=8, frameon=False, markerscale=3, loc='upper right')
ax.set_title('MMR Status', fontsize=10, fontweight='bold')
ax.set_xticks([]); ax.set_yticks([])
ax.spines[['top','right','bottom','left']].set_visible(False)

# 3. UMAP — attention score
ax = fig.add_subplot(gs[0, 2])
vmax_att = np.percentile(att_scores, 98)
sc3 = ax.scatter(xy[:, 0], xy[:, 1], c=att_scores, cmap='YlOrRd',
                  vmin=0, vmax=vmax_att, s=1.2, alpha=0.6,
                  linewidths=0, rasterized=True)
plt.colorbar(sc3, ax=ax, shrink=0.7, pad=0.02, label='Attention score')
ax.set_title('ABMIL Attention Score\n(Immune-only bags)', fontsize=10, fontweight='bold')
ax.set_xticks([]); ax.set_yticks([])
ax.spines[['top','right','bottom','left']].set_visible(False)

# 4&5. UMAP — attention per MMR group
for col_idx, (mmr_grp, cmap_) in enumerate([('MMRp','YlOrRd'), ('MMRd','YlGnBu')]):
    ax = fig.add_subplot(gs[1, col_idx])
    mask_fg = adata.obs['MMRStatus'] == mmr_grp
    mask_bg = ~mask_fg
    ax.scatter(xy[mask_bg, 0], xy[mask_bg, 1], c='#e0e0e0',
               s=0.8, alpha=0.25, linewidths=0, rasterized=True)
    sc_ = ax.scatter(xy[mask_fg, 0], xy[mask_fg, 1],
                      c=att_scores[mask_fg], cmap=cmap_,
                      vmin=0, vmax=vmax_att,
                      s=1.5, alpha=0.75, linewidths=0, rasterized=True)
    plt.colorbar(sc_, ax=ax, shrink=0.7, pad=0.02,
                  label=f'Attention ({mmr_grp})')
    ax.set_title(f'Attention — {mmr_grp} cells', fontsize=10, fontweight='bold')
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[['top','right','bottom','left']].set_visible(False)

# 6. CV AUC bar
ax = fig.add_subplot(gs[1, 2])
fold_labels = [f'F{i+1}' for i in range(N_FOLDS)]
bar_colors  = ['#E64B35' if a == max(fold_aucs) else '#8491B4' for a in fold_aucs]
bars = ax.bar(fold_labels, fold_aucs, color=bar_colors, alpha=0.85,
               edgecolor='white', width=0.6)
ax.axhline(np.mean(fold_aucs), color='black', lw=1.5, ls='--',
           label=f'Mean={np.mean(fold_aucs):.3f}±{np.std(fold_aucs):.3f}')
ax.axhline(0.5, color='grey', lw=0.8, ls=':', alpha=0.6, label='Chance')
for b, a in zip(bars, fold_aucs):
    ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.01,
            f'{a:.3f}', ha='center', va='bottom', fontsize=8)
ax.set_ylim(0.3, 1.1)
ax.set_ylabel('ROC-AUC', fontsize=9)
ax.set_title('5-Fold CV (Immune-only bags)\nMMRp vs MMRd',
              fontsize=10, fontweight='bold')
ax.legend(fontsize=7.5, frameon=False)
ax.spines[['top', 'right']].set_visible(False)

# 7. Violin — attention by subtype (MMRp bags, sorted by mean)
ax = fig.add_subplot(gs[2, :2])

subtype_order = stats_mmrp['subtype'].tolist()
# 仅保留样本量足够的亚群
violin_data = [
    mmrp_obs[mmrp_obs['clMidwayPr'] == s]['attention_score'].values
    for s in subtype_order
]
positions = list(range(len(subtype_order)))

parts = ax.violinplot(violin_data, positions=positions, widths=0.7,
                       showmedians=True, showextrema=False)
for i, (body, st) in enumerate(zip(parts['bodies'], subtype_order)):
    col = SUBTYPE_COLORS.get(st, '#aaaaaa')
    body.set_facecolor(col)
    body.set_alpha(0.65)
parts['cmedians'].set_color('black')
parts['cmedians'].set_linewidth(1.5)

# 高亮 Pre-driver 区域
for st in ['Macro', 'TCD8']:
    if st in subtype_order:
        xi = subtype_order.index(st)
        ax.axvspan(xi - 0.45, xi + 0.45, alpha=0.08, color='#B22222', zorder=0)

# FDR 星号标注
if len(pairs_df) > 0:
    y_ann = max(v.max() for v in violin_data if len(v) > 0) * 1.03
    for _, row in pairs_df.iterrows():
        st = row['subtype']
        if st in subtype_order:
            xi = subtype_order.index(st)
            fdr = row['fdr']
            sig = '***' if fdr < 0.001 else '**' if fdr < 0.01 else '*' if fdr < 0.05 else 'ns'
            if sig != 'ns':
                ax.text(xi, y_ann, sig, ha='center', fontsize=9,
                        fontweight='bold', color=SUBTYPE_COLORS.get(st,'black'))

# 均值点
for i, data in enumerate(violin_data):
    if len(data) > 0:
        ax.scatter([i], [np.mean(data)], marker='D', color='white',
                    edgecolors='black', s=30, zorder=5, linewidths=0.8)

ax.set_xticks(positions)
ax.set_xticklabels(subtype_order, rotation=30, ha='right', fontsize=9)
ax.set_ylabel('Attention Score', fontsize=9)
ax.set_title('Per-cell Attention Score by Immune Subtype — MMRp patients\n'
             '(sorted by mean attention; ◆ = mean; shaded = Pre-driver subtypes)',
              fontsize=9.5, fontweight='bold')
ax.spines[['top', 'right']].set_visible(False)

# 8. MMRp vs MMRd heatmap — mean attention by subtype
ax = fig.add_subplot(gs[2, 2])

all_subtypes_h = sorted(
    set(stats_mmrp['subtype'].tolist() + stats_mmrd['subtype'].tolist())
)
# Build matrix
heat_data = np.full((len(all_subtypes_h), 2), np.nan)
mmrp_dict = dict(zip(stats_mmrp['subtype'], stats_mmrp['mean_att']))
mmrd_dict = dict(zip(stats_mmrd['subtype'], stats_mmrd['mean_att']))
for i, st in enumerate(all_subtypes_h):
    heat_data[i, 0] = mmrp_dict.get(st, np.nan)
    heat_data[i, 1] = mmrd_dict.get(st, np.nan)

# normalize per-column for display (z-score within MMR group)
heat_z = heat_data.copy()
for col in range(2):
    col_vals = heat_data[:, col]
    valid = ~np.isnan(col_vals)
    heat_z[valid, col] = (col_vals[valid] - col_vals[valid].mean()) / \
                          (col_vals[valid].std() + 1e-10)
heat_z = np.clip(heat_z, -2.5, 2.5)

im = ax.imshow(heat_z.T, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5)
plt.colorbar(im, ax=ax, shrink=0.55, pad=0.04, label='z-score')
ax.set_xticks(range(len(all_subtypes_h)))
ax.set_xticklabels(all_subtypes_h, rotation=45, ha='right', fontsize=7.8)
ax.set_yticks([0, 1])
ax.set_yticklabels(['MMRp', 'MMRd'], fontsize=9)
ax.set_title('Attention z-score\n(MMRp vs MMRd)', fontsize=9, fontweight='bold')
ax.spines[['top','right','bottom','left']].set_visible(False)

# 高亮 Macro & TCD8
for st in ['Macro', 'TCD8']:
    if st in all_subtypes_h:
        xi = all_subtypes_h.index(st)
        ax.axvline(xi, color='#B22222', lw=1.5, alpha=0.5, zorder=5)
        ax.get_xticklabels()[xi].set_color('#B22222')
        ax.get_xticklabels()[xi].set_fontweight('bold')

# Global title
mean_auc = np.mean(fold_aucs)
std_auc  = np.std(fold_aucs)
fig.text(0.5, 0.975,
         f'Figure 3B — Immune-only ABMIL: Attention Analysis  '
         f'(5-fold CV AUC={mean_auc:.3f}±{std_auc:.3f})',
         fontsize=13, fontweight='bold', ha='center', va='top')
fig.text(0.5, 0.965,
         'Bags restricted to TNKILC + Myeloid cells  |  '
         'Shaded subtypes: Phase 2 Pre-driver candidates  |  '
         '★ FDR significance vs rest',
         fontsize=8.5, ha='center', va='top', color='grey', style='italic')

out_name = 'figure3B_mil_immune_attention'
fig.savefig(f'{FIG_DIR}/{out_name}.png', dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(f'{FIG_DIR}/{out_name}.pdf',           bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'  Saved: {FIG_DIR}/{out_name}.png/.pdf')

# 保存 attention scores 回 full data
adata.write_h5ad('processed_data/GSE178341_immune_mil.h5ad', compression='gzip')
att_out = adata.obs[['PID', 'MMRStatus', 'clTopLevel', 'clMidwayPr',
                      'cl295v11SubShort', 'attention_score']].copy()
att_out.to_csv(f'{RES_DIR}/immune_attention_scores.csv')

print('\n' + '='*60)
print(f'DONE — Immune-only MIL')
print(f'  CV AUC: {mean_auc:.4f} ± {std_auc:.4f}')
if 'Macro' in tests: print(f'  Macro vs others (MMRp): p={tests["Macro"]:.4e}')
if 'TCD8'  in tests: print(f'  TCD8  vs others (MMRp): p={tests["TCD8"]:.4e}')
