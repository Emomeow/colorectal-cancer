#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 3 Panel E & F
Panel E: Augur Predictability — cell type sensitivity to MMR status
Panel F: Transfer Validation — cold/hot scoring + tumor vs paired normal attention
"""

import scanpy as sc
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.decomposition import PCA
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import scipy.sparse as sp
import io, sys, os, warnings
warnings.filterwarnings('ignore')

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

torch.manual_seed(42)
np.random.seed(42)

plt.rcParams.update({
    'font.family':  'Arial',
    'pdf.fonttype': 42,
    'ps.fonttype':  42,
    'axes.spines.top':   False,
    'axes.spines.right': False,
})

FIG_DIR = 'figures/phase3'
RES_DIR = 'phase3_results/transfer'
os.makedirs(FIG_DIR, exist_ok=True)

PRE_DRIVER = {'Macro', 'TCD8'}
SUBTYPE_COLORS = {
    'Macro': '#E64B35', 'TCD8': '#3C5488', 'TCD4': '#00A087',
    'Treg': '#F39B7F',  'NK':   '#8491B4', 'DC':   '#91D1C2',
    'Tgd':  '#B09C85',  'Other':'#DDDDDD',
}
GSE132_SUBTYPE_MAP = {
    'CD8+ T cells': 'TCD8', 'CD4+ T cells': 'TCD4',
    'T helper 17 cells': 'TCD4', 'T follicular helper cells': 'TCD4',
    'Regulatory T cells': 'Treg', 'gamma delta T cells': 'Tgd',
    'NK cells': 'NK', 'SPP1+': 'Macro', 'Pro-inflammatory': 'Macro',
    'cDC': 'DC', 'Proliferating': 'Other', 'Unknown': 'Other',
}
N_COMPS   = 50
MAX_CELLS = 1500
HIDDEN_DIM = 256
ATT_DIM    = 128
DROPOUT    = 0.25

# ══════════════════════════════════════════════════════════════════════════
#  Panel E: Augur Predictability (load pre-computed results)
# ══════════════════════════════════════════════════════════════════════════
print('Plotting Panel E (Augur) ...')
augur_df = pd.read_csv('phase3_results/augur/augur_auc.csv')
augur_df = augur_df.sort_values('auc', ascending=True).reset_index(drop=True)

fig_e, ax = plt.subplots(figsize=(6.5, max(4, len(augur_df) * 0.55 + 1.2)))
y_pos = np.arange(len(augur_df))

bar_colors = [SUBTYPE_COLORS.get(ct, '#CCCCCC') if ct in PRE_DRIVER else '#CCCCCC'
              for ct in augur_df['cell_type']]

ax.axvline(0.5, color='black', lw=1.0, ls='--', alpha=0.5)
ax.axvspan(0.5, 1.0, alpha=0.04, color='#F39B7F')

for i, row in augur_df.iterrows():
    ax.plot([0.5, row['auc']], [i, i], color=bar_colors[i], lw=1.8, zorder=2)
ax.scatter(augur_df['auc'], y_pos,
           s=[120 if ct in PRE_DRIVER else 70 for ct in augur_df['cell_type']],
           c=bar_colors, zorder=3, edgecolors='white', linewidths=0.5)
ax.errorbar(augur_df['auc'], y_pos, xerr=augur_df['auc_std'],
            fmt='none', color='gray', elinewidth=0.8, capsize=2, zorder=2)

for i, row in augur_df.iterrows():
    ax.text(row['auc'] + row['auc_std'] + 0.008, i,
            f'{row["auc"]:.3f}', va='center', ha='left', fontsize=8,
            color='black' if row['cell_type'] in PRE_DRIVER else '#666666')

ax.set_yticks(y_pos)
ax.set_yticklabels(augur_df['cell_type'], fontsize=10)
for tick, ct in zip(ax.get_yticklabels(), augur_df['cell_type']):
    if ct in PRE_DRIVER:
        tick.set_color(SUBTYPE_COLORS.get(ct, '#E64B35'))
        tick.set_fontweight('bold')
for i, row in augur_df.iterrows():
    if row['cell_type'] in PRE_DRIVER:
        ax.text(1.02, i, 'Pre-driver', transform=ax.get_yaxis_transform(),
                va='center', ha='left', fontsize=7.5,
                color=SUBTYPE_COLORS.get(row['cell_type'], '#E64B35'),
                fontstyle='italic')

ax.set_xlabel('Augur AUC  (MMRp vs MMRd)', fontsize=11)
ax.set_title('Cell Type Predictability (Augur)\nImmune Microenvironment · GSE178341',
             fontsize=12, fontweight='bold', pad=12)
ax.set_xlim(0.35, 1.05)
ax.set_clip_on(False)
patches = [
    mpatches.Patch(color='#E64B35', label='Pre-driver: Macro (TAM)'),
    mpatches.Patch(color='#3C5488', label='Pre-driver: TCD8 (Tex)'),
    mpatches.Patch(color='#CCCCCC', label='Other immune subtypes'),
    plt.Line2D([0], [0], ls='--', color='black', alpha=0.5, label='Random baseline'),
]
ax.legend(handles=patches, fontsize=8, loc='lower right', framealpha=0.7, edgecolor='lightgray')
plt.tight_layout()
fig_e.savefig(f'{FIG_DIR}/figure3E_augur_predictability.png', dpi=200, bbox_inches='tight')
fig_e.savefig(f'{FIG_DIR}/figure3E_augur_predictability.pdf', bbox_inches='tight')
plt.close(fig_e)
print('  Saved figure3E_augur_predictability.png/.pdf')

# ══════════════════════════════════════════════════════════════════════════
#  Prepare model + PCA for Panel F inference
# ══════════════════════════════════════════════════════════════════════════
print('\nPreparing model for Panel F ...')

# Load GSE178341, fit PCA on common HVGs
adata178 = sc.read_h5ad('processed_data/GSE178341_immune_mil.h5ad')
merged_genes = set(sc.read_h5ad('processed_data/merged_annotated.h5ad', backed='r').var_names)
hvg_all  = adata178.var_names[adata178.var['highly_variable']].tolist()
common_hvg = [g for g in hvg_all if g in merged_genes]
print(f'  Common HVGs: {len(common_hvg)}')

X178 = adata178[:, common_hvg].X
if sp.issparse(X178): X178 = X178.toarray()
X178 = X178.astype(np.float32)

pca_model = PCA(n_components=N_COMPS, svd_solver='arpack', random_state=42)
pca_model.fit(X178)
print(f'  PCA variance explained: {pca_model.explained_variance_ratio_.sum():.3f}')

# ABMIL model definition
class GatedABMIL(nn.Module):
    def __init__(self):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(N_COMPS, HIDDEN_DIM), nn.ReLU(), nn.Dropout(DROPOUT),
            nn.Linear(HIDDEN_DIM, HIDDEN_DIM), nn.ReLU(),
        )
        self.att_V = nn.Sequential(nn.Linear(HIDDEN_DIM, ATT_DIM), nn.Tanh())
        self.att_U = nn.Sequential(nn.Linear(HIDDEN_DIM, ATT_DIM), nn.Sigmoid())
        self.att_w = nn.Linear(ATT_DIM, 1, bias=False)
        self.classifier = nn.Sequential(nn.Dropout(DROPOUT), nn.Linear(HIDDEN_DIM, 1))

    def forward(self, x, return_attention=False):
        H = self.encoder(x)
        A = F.softmax(self.att_w(self.att_V(H) * self.att_U(H)), dim=0)
        M = (A * H).sum(dim=0, keepdim=True)
        logit = self.classifier(M).squeeze()
        return (logit, A.squeeze()) if return_attention else logit

model = GatedABMIL()
model.load_state_dict(torch.load(f'{RES_DIR}/abmil_transfer.pt', weights_only=True))
model.eval()
print('  Model loaded.')

def build_feature_matrix(adata_sub, common_hvg, pca_model):
    """Normalize counts → log1p → PCA projection."""
    # Check if already normalized
    vals = adata_sub.X[:5, :].toarray() if sp.issparse(adata_sub.X) else adata_sub.X[:5, :]
    if vals.max() > 100:  # raw counts
        sc.pp.normalize_total(adata_sub, target_sum=1e4)
        sc.pp.log1p(adata_sub)
    X = np.zeros((adata_sub.n_obs, len(common_hvg)), dtype=np.float32)
    gene_idx = {g: i for i, g in enumerate(adata_sub.var_names)}
    for j, gene in enumerate(common_hvg):
        if gene in gene_idx:
            col = adata_sub.X[:, gene_idx[gene]]
            if sp.issparse(col): col = col.toarray().ravel()
            else: col = np.asarray(col).ravel()
            X[:, j] = col
    return pca_model.transform(X).astype(np.float32)

@torch.no_grad()
def infer_attention(model, X_pca, pids, bags):
    att = np.zeros(len(X_pca), dtype=np.float32)
    probs = {}
    for pid in np.unique(pids):
        idx = np.where(pids == pid)[0]
        x   = torch.tensor(X_pca[idx])
        logit, A = model(x, return_attention=True)
        att[idx] = A.numpy()
        probs[pid] = torch.sigmoid(logit).item()
    return att, probs

# ══════════════════════════════════════════════════════════════════════════
#  Load GSE132465: both Tumor + Normal immune cells
# ══════════════════════════════════════════════════════════════════════════
print('\nLoading GSE132465 ...')
adata_merged = sc.read_h5ad('processed_data/merged_annotated.h5ad')
IMMUNE_TYPES = ['T cells', 'Myeloids', 'NK cells']

gse = adata_merged[(adata_merged.obs['batch'] == 'GSE132465') &
                    adata_merged.obs['unified_cell_type'].isin(IMMUNE_TYPES)].copy()
del adata_merged

# Split tumor vs normal
tumor  = gse[gse.obs['Class'] == 'Tumor'].copy()
normal = gse[gse.obs['Class'] == 'Normal'].copy()
print(f'  Tumor immune: {tumor.n_obs:,}  Normal immune: {normal.n_obs:,}')

PAIRED_PATIENTS = sorted(set(tumor.obs['Patient'].unique()) &
                           set(normal.obs['Patient'].unique()))
print(f'  Paired patients (tumor+normal): {len(PAIRED_PATIENTS)}')

# Project + infer
print('  Processing tumor cells ...')
X_tumor_pca  = build_feature_matrix(tumor, common_hvg, pca_model)
att_tumor, probs_tumor = infer_attention(model, X_tumor_pca,
                                          tumor.obs['Patient'].values, {})
tumor.obs['attention'] = att_tumor
tumor.obs['subtype']   = tumor.obs['unified_cell_subtype'].map(GSE132_SUBTYPE_MAP).fillna('Other')

print('  Processing normal cells ...')
X_normal_pca = build_feature_matrix(normal, common_hvg, pca_model)
att_normal, probs_normal = infer_attention(model, X_normal_pca,
                                             normal.obs['Patient'].values, {})
normal.obs['attention'] = att_normal
normal.obs['subtype']   = normal.obs['unified_cell_subtype'].map(GSE132_SUBTYPE_MAP).fillna('Other')

# Cold/hot scores (tumor)
pred_df = pd.read_csv(f'{RES_DIR}/gse132465_predictions.csv')
print(f'\n  Cold/hot prediction range: {pred_df["prob_MMRd"].min():.3f} – {pred_df["prob_MMRd"].max():.3f}')

# ── Signature scoring：肿瘤特异性验证（不依赖 ABMIL out-of-distribution 行为）──
print('\n  Computing signature scores for tumor vs normal ...')

# 肿瘤免疫抑制 TAM 相关（merged data 中仅有以下基因可用）
TAM_PROTUMOR = ['VEGFA', 'CD163', 'HMOX1', 'LGALS3', 'STAT3', 'CD68', 'CD14']
T_EXHAUSTION  = ['PDCD1', 'LAG3', 'TIGIT', 'HAVCR2', 'CTLA4',
                  'TOX', 'ENTPD1', 'BATF']

def score_sig(adata_sub, genes, name):
    avail = [g for g in genes if g in adata_sub.var_names]
    if len(avail) >= 3:
        sc.tl.score_genes(adata_sub, avail, score_name=name, use_raw=False)
    else:
        adata_sub.obs[name] = 0.0

score_sig(tumor,  TAM_PROTUMOR, 'TAM_ProTumor')
score_sig(tumor,  T_EXHAUSTION,  'T_Exhaustion')
score_sig(normal, TAM_PROTUMOR, 'TAM_ProTumor')
score_sig(normal, T_EXHAUSTION,  'T_Exhaustion')

# 配对患者 per-patient mean signature score
FOCUS_CONFIGS = [
    ('TCD8',  'T_Exhaustion'),
    ('TCD4',  'T_Exhaustion'),   # control
]
paired_records = []
for pid in PAIRED_PATIENTS:
    for sub, sig in FOCUS_CONFIGS:
        t_vals = tumor.obs.loc[(tumor.obs['Patient'] == pid) &
                                (tumor.obs['subtype'] == sub), sig].dropna().values
        n_vals = normal.obs.loc[(normal.obs['Patient'] == pid) &
                                 (normal.obs['subtype'] == sub), sig].dropna().values
        if len(t_vals) >= 5 and len(n_vals) >= 3:
            paired_records.append({
                'patient': pid, 'subtype': sub, 'signature': sig,
                'tumor_mean': t_vals.mean(), 'normal_mean': n_vals.mean(),
                'n_tumor': len(t_vals), 'n_normal': len(n_vals),
            })

paired_df = pd.DataFrame(paired_records)
paired_df.to_csv(f'{RES_DIR}/paired_tumor_normal_sig.csv', index=False)
print('  Signature scores (paired tumor vs normal):')
print(paired_df.groupby(['subtype', 'signature'])[['tumor_mean', 'normal_mean']].mean().round(3).to_string())

# Wilcoxon signed-rank tests
print('\n  Wilcoxon signed-rank tests (tumor > normal):')
stat_results = []
for sub, sig in FOCUS_CONFIGS:
    sub_df = paired_df[(paired_df['subtype'] == sub) & (paired_df['signature'] == sig)]
    if len(sub_df) >= 3:
        stat, p = wilcoxon(sub_df['tumor_mean'], sub_df['normal_mean'], alternative='greater')
        print(f'  {sub:8s} {sig:20s}: n={len(sub_df)}, p={p:.4e}')
        stat_results.append({'subtype': sub, 'signature': sig, 'n_pairs': len(sub_df), 'pval': p})

# ══════════════════════════════════════════════════════════════════════════
#  Panel F: (left) Cold/Hot scoring  |  (right) Tumor vs Normal attention
# ══════════════════════════════════════════════════════════════════════════
print('\nPlotting Panel F ...')

fig_f = plt.figure(figsize=(14.5, 5.5))
gs = gridspec.GridSpec(1, 3, width_ratios=[1.3, 0.85, 0.85], wspace=0.38,
                       left=0.07, right=0.97, top=0.87, bottom=0.14)

# ─── Left: Cold / Hot Scoring Bar Chart ─────────────────────────────────
ax_l = fig_f.add_subplot(gs[0])

pred_sorted = pred_df.sort_values('prob_MMRd').reset_index(drop=True)
probs = pred_sorted['prob_MMRd'].values
patients = pred_sorted['patient'].values
y = np.arange(len(pred_sorted))

# Colormap: cold (blue) → hot (red) at threshold 0.5
cmap = matplotlib.colormaps.get_cmap('RdYlBu_r')
bar_colors = [cmap(p) for p in probs]

bars = ax_l.barh(y, probs, color=bar_colors, height=0.65,
                  edgecolor='white', linewidth=0.3)
ax_l.axvline(0.5, color='black', lw=1.2, ls='--', alpha=0.7,
             label='Classification threshold (0.5)')

# Threshold labels
ax_l.text(0.22, len(pred_sorted), 'Cold (MMRp-like)', ha='center', va='bottom',
          fontsize=8.5, color='#4DBBD5', fontstyle='italic')
ax_l.text(0.78, len(pred_sorted), 'Hot (MMRd-like)', ha='center', va='bottom',
          fontsize=8.5, color='#E64B35', fontstyle='italic')

# Score labels
for i, (p, pid) in enumerate(zip(probs, patients)):
    ax_l.text(max(p + 0.02, 0.04), i, f'{p:.2f}', va='center', ha='left',
              fontsize=6.8, color='#333333')

ax_l.set_yticks(y)
ax_l.set_yticklabels(patients, fontsize=8.5)
ax_l.set_xlabel('Predicted MMRd Probability ("热" Score)', fontsize=10)
ax_l.set_xlim(0, 1.15)
ax_l.set_title('Cold / Hot Tumor Scoring\nGSE132465 · 23 CRC Patients',
               fontsize=11, fontweight='bold')
ax_l.legend(fontsize=8, loc='lower right', framealpha=0.8)

# Colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
sm.set_array([])
cb = plt.colorbar(sm, ax=ax_l, shrink=0.45, aspect=15, pad=0.02, location='right')
cb.set_label('Prob. MMRd', fontsize=8)
cb.ax.tick_params(labelsize=7)

# ─── Right: Tumor vs Normal — Signature Score Paired Plot ────────────────
ax_r = fig_f.add_subplot(gs[1])

# 3 panels: Macro/TAM_ProTumor, TCD8/T_Exhaustion, TCD4/T_Exhaustion (control)
plot_configs = [
    ('TCD8',  'T_Exhaustion',  'TCD8\nT-Exhaustion\n(Pre-driver)'),
    ('TCD4',  'T_Exhaustion',  'TCD4\nT-Exhaustion\n(control)'),
]
n_sub  = len(plot_configs)
x_base = np.arange(n_sub)
W      = 0.25

for xi, (sub, sig, label) in enumerate(plot_configs):
    sub_df = paired_df[(paired_df['subtype'] == sub) & (paired_df['signature'] == sig)]
    if len(sub_df) == 0:
        continue
    t_vals = sub_df['tumor_mean'].values
    n_vals = sub_df['normal_mean'].values

    c_tumor  = SUBTYPE_COLORS.get(sub, '#888888')
    c_normal = '#BBBBBB'

    # Paired connecting lines
    for tv, nv in zip(t_vals, n_vals):
        ax_r.plot([xi - W * 0.5, xi + W * 0.5], [nv, tv],
                  color='#AAAAAA', lw=0.8, alpha=0.55, zorder=1)

    # Dots: Normal (left) and Tumor (right)
    rng = np.random.RandomState(xi + 7)
    jitter = rng.uniform(-0.045, 0.045, len(t_vals))
    ax_r.scatter(xi - W * 0.5 + jitter, n_vals, s=30, color=c_normal,
                 zorder=3, edgecolors='gray', linewidths=0.5, alpha=0.9)
    ax_r.scatter(xi + W * 0.5 + jitter, t_vals, s=38, color=c_tumor,
                 zorder=3, edgecolors='white', linewidths=0.5, alpha=0.92)

    # Mean bars
    ax_r.plot([xi - W * 0.5 - 0.08, xi - W * 0.5 + 0.08],
              [n_vals.mean()] * 2, color='#666666', lw=2.5, zorder=5)
    ax_r.plot([xi + W * 0.5 - 0.08, xi + W * 0.5 + 0.08],
              [t_vals.mean()] * 2, color=c_tumor, lw=2.5, zorder=5)

    # Statistical annotation
    sr = next((r for r in stat_results
               if r['subtype'] == sub and r['signature'] == sig), None)
    if sr:
        p = sr['pval']
        stars = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'ns'))
        y_top = max(t_vals.max(), n_vals.max())
        y_ann = y_top + abs(y_top) * 0.20 + 0.05
        ax_r.plot([xi - W * 0.5, xi - W * 0.5, xi + W * 0.5, xi + W * 0.5],
                  [y_ann * 0.93, y_ann, y_ann, y_ann * 0.93],
                  lw=1.0, color='black')
        color_s = '#E64B35' if stars != 'ns' else 'gray'
        ax_r.text(xi, y_ann * 1.01, stars, ha='center', va='bottom',
                  fontsize=11, color=color_s, fontweight='bold')

# N/T labels below x-axis
for xi in range(n_sub):
    ax_r.text(xi - W * 0.5, -0.08, 'N', ha='center', va='top',
              fontsize=8, color='gray', transform=ax_r.get_xaxis_transform())
    ax_r.text(xi + W * 0.5, -0.08, 'T', ha='center', va='top',
              fontsize=8, color='black', transform=ax_r.get_xaxis_transform())

ax_r.set_xticks(x_base)
ax_r.set_xticklabels([cfg[2] for cfg in plot_configs], fontsize=9)
ax_r.set_ylabel('Mean Exhaustion Score', fontsize=10)
ax_r.set_title('T Cell Exhaustion:\nTumor vs Normal · n=10 paired',
               fontsize=10.5, fontweight='bold')
ax_r.set_xlim(-0.55, n_sub - 0.45)

leg_handles = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#BBBBBB',
               markersize=8, markeredgecolor='gray', label='Normal (N)'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#3C5488',
               markersize=8, label='Tumor (T)'),
]
ax_r.legend(handles=leg_handles, fontsize=8, loc='upper right',
            framealpha=0.85, edgecolor='lightgray')

# ─── Third panel: Relative attention (Pre-driver / all) in Tumor vs Normal ─
ax_r2 = fig_f.add_subplot(gs[2])

# Per-patient: mean_att(Pre-driver) / mean_att(all) in Tumor vs Normal
rel_records = []
for pid in PAIRED_PATIENTS:
    for ct in ['Macro', 'TCD8']:
        t_pd  = tumor.obs.loc[(tumor.obs['Patient'] == pid) &
                               (tumor.obs['subtype'] == ct), 'attention'].values
        t_all = tumor.obs.loc[tumor.obs['Patient'] == pid, 'attention'].values
        n_pd  = normal.obs.loc[(normal.obs['Patient'] == pid) &
                                (normal.obs['subtype'] == ct), 'attention'].values
        n_all = normal.obs.loc[normal.obs['Patient'] == pid, 'attention'].values
        if len(t_pd) < 5 or len(n_pd) < 3:
            continue
        t_ratio = t_pd.mean() / (t_all.mean() + 1e-10)
        n_ratio = n_pd.mean() / (n_all.mean() + 1e-10)
        rel_records.append({'patient': pid, 'subtype': ct,
                             'tumor_ratio': t_ratio, 'normal_ratio': n_ratio})

rel_df = pd.DataFrame(rel_records)
print('\n  Relative attention (Pre-driver / all):')
if len(rel_df):
    print(rel_df.groupby('subtype')[['tumor_ratio', 'normal_ratio']].mean().round(3).to_string())
    rel_df.to_csv(f'{RES_DIR}/relative_attention_paired.csv', index=False)

# Plot side-by-side for Macro and TCD8
rel_order = ['TCD8']
rel_labels = ['TCD8\n(Pre-driver)']
x_rel = np.arange(len(rel_order))
W_r = 0.25

for xi, (sub, label) in enumerate(zip(rel_order, rel_labels)):
    sub_df = rel_df[rel_df['subtype'] == sub] if len(rel_df) else pd.DataFrame()
    if len(sub_df) == 0:
        continue
    t_vals = sub_df['tumor_ratio'].values
    n_vals = sub_df['normal_ratio'].values

    c_sub = SUBTYPE_COLORS.get(sub, '#888888')
    for tv, nv in zip(t_vals, n_vals):
        ax_r2.plot([xi - W_r * 0.5, xi + W_r * 0.5], [nv, tv],
                   color='#AAAAAA', lw=0.8, alpha=0.55, zorder=1)

    rng = np.random.RandomState(xi + 11)
    jit = rng.uniform(-0.04, 0.04, len(t_vals))
    ax_r2.scatter(xi - W_r * 0.5 + jit, n_vals, s=30, color='#BBBBBB',
                  zorder=3, edgecolors='gray', linewidths=0.5)
    ax_r2.scatter(xi + W_r * 0.5 + jit, t_vals, s=38, color=c_sub,
                  zorder=3, edgecolors='white', linewidths=0.5)
    ax_r2.plot([xi - W_r * 0.5 - 0.07, xi - W_r * 0.5 + 0.07],
               [n_vals.mean()] * 2, color='#666666', lw=2.5, zorder=5)
    ax_r2.plot([xi + W_r * 0.5 - 0.07, xi + W_r * 0.5 + 0.07],
               [t_vals.mean()] * 2, color=c_sub, lw=2.5, zorder=5)

    if len(t_vals) >= 3:
        try:
            _, p_r = wilcoxon(t_vals, n_vals, alternative='two-sided')
            stars_r = '**' if p_r < 0.01 else ('*' if p_r < 0.05 else 'ns')
            y_top = max(t_vals.max(), n_vals.max())
            y_ann = y_top * 1.18 + 0.05
            ax_r2.plot([xi - W_r * 0.5, xi - W_r * 0.5, xi + W_r * 0.5, xi + W_r * 0.5],
                       [y_ann * 0.94, y_ann, y_ann, y_ann * 0.94], lw=1.0, color='black')
            ax_r2.text(xi, y_ann * 1.02, stars_r, ha='center', va='bottom',
                       fontsize=10, color='#E64B35' if stars_r != 'ns' else 'gray',
                       fontweight='bold')
        except:
            pass

ax_r2.axhline(1.0, color='black', lw=0.8, ls=':', alpha=0.5, label='Equal to mean')
ax_r2.set_xticks(x_rel)
ax_r2.set_xticklabels(rel_labels, fontsize=9)
ax_r2.set_ylabel('Relative Attention\n(Pre-driver / All immune)', fontsize=9.5)
ax_r2.set_title('TCD8 Relative Attention:\nTumor vs Paired Normal (n=10)',
                fontsize=10.5, fontweight='bold')
ax_r2.set_xlim(-0.6, 0.6)
for xi in range(len(rel_order)):
    ax_r2.text(xi - W_r * 0.5, -0.06, 'N', ha='center', va='top', fontsize=8,
               color='gray', transform=ax_r2.get_xaxis_transform())
    ax_r2.text(xi + W_r * 0.5, -0.06, 'T', ha='center', va='top', fontsize=8,
               color='black', transform=ax_r2.get_xaxis_transform())

fig_f.suptitle('Transfer Validation — Cold/Hot Classification & Pre-driver Tumor Specificity',
               fontsize=12, fontweight='bold')
fig_f.savefig(f'{FIG_DIR}/figure3F_transfer_validation.png', dpi=200, bbox_inches='tight')
fig_f.savefig(f'{FIG_DIR}/figure3F_transfer_validation.pdf', bbox_inches='tight')
plt.close(fig_f)
print('  Saved figure3F_transfer_validation.png/.pdf')

print('\nDone.')
