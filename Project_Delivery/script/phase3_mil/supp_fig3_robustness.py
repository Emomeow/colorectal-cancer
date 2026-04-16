#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplementary Figure 3 — ABMIL Hyperparameter Robustness
GSE178341 immune-only bags · 5-fold CV OOF AUC across parameter sweeps
"""

import scanpy as sc
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import Adam
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
})

FIG_DIR = 'figures/phase3'
os.makedirs(FIG_DIR, exist_ok=True)

# ── Baseline hyperparameters ────────────────────────────────────────────────
BASELINE = dict(
    N_PCA      = 50,
    HIDDEN_DIM = 256,
    ATT_DIM    = 128,
    MAX_CELLS  = 1500,
    DROPOUT    = 0.25,
    LR         = 1e-4,
)
N_FOLDS   = 5
EPOCHS    = 80       # slightly reduced for sweep speed
PATIENCE  = 15
WD        = 1e-5

# ── Sweep grid (one param at a time, others = baseline) ────────────────────
SWEEPS = {
    'N_PCA':      [10, 20, 30, 40, 50],
    'HIDDEN_DIM': [64, 128, 256, 512],
    'ATT_DIM':    [32, 64, 128, 256],
    'MAX_CELLS':  [500, 750, 1000, 1500, 2000],
    'DROPOUT':    [0.0, 0.10, 0.25, 0.40, 0.50],
    'LR':         [1e-5, 5e-5, 1e-4, 5e-4, 1e-3],
}

PARAM_LABELS = {
    'N_PCA':      'N PCA Components',
    'HIDDEN_DIM': 'Hidden Dim',
    'ATT_DIM':    'Attention Dim',
    'MAX_CELLS':  'Max Cells per Bag',
    'DROPOUT':    'Dropout Rate',
    'LR':         'Learning Rate',
}

# ══════════════════════════════════════════════════════════════════════════
#  Load data
# ══════════════════════════════════════════════════════════════════════════
print('Loading data ...')
adata = sc.read_h5ad('processed_data/GSE178341_immune_mil.h5ad')
print(f'  {adata.n_obs:,} cells, {adata.obs["PID"].nunique()} patients')

# Full PCA matrix (100 dims — we'll slice as needed)
pca_full = adata.obsm['X_pca'].astype('float32')   # use all stored dims
print(f'  PCA dims available: {pca_full.shape[1]}')
pids     = adata.obs['PID'].values
mmr_lbl  = (adata.obs['MMRStatus'] == 'MMRd').astype(int).values

unique_pids = adata.obs['PID'].unique()
patient_meta = {}
for pid in unique_pids:
    idx = np.where(pids == pid)[0]
    patient_meta[pid] = {'idx': idx, 'label': int(mmr_lbl[idx[0]])}

all_pids   = np.array(list(patient_meta.keys()))
all_labels = np.array([patient_meta[p]['label'] for p in all_pids])
print(f'  Bags: {len(all_pids)}  (MMRd={all_labels.sum()}, MMRp={(1-all_labels).sum()})')

skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)
folds = list(skf.split(all_pids, all_labels))


# ── Model factory ──────────────────────────────────────────────────────────
def build_model(n_pca, hidden_dim, att_dim, dropout):
    class GatedABMIL(nn.Module):
        def __init__(self):
            super().__init__()
            self.encoder = nn.Sequential(
                nn.Linear(n_pca, hidden_dim), nn.ReLU(), nn.Dropout(dropout),
                nn.Linear(hidden_dim, hidden_dim), nn.ReLU(),
            )
            self.att_V = nn.Sequential(nn.Linear(hidden_dim, att_dim), nn.Tanh())
            self.att_U = nn.Sequential(nn.Linear(hidden_dim, att_dim), nn.Sigmoid())
            self.att_w = nn.Linear(att_dim, 1, bias=False)
            self.classifier = nn.Sequential(nn.Dropout(dropout),
                                            nn.Linear(hidden_dim, 1))
        def forward(self, x):
            H = self.encoder(x)
            A = F.softmax(self.att_w(self.att_V(H) * self.att_U(H)), dim=0)
            return self.classifier((A * H).sum(dim=0, keepdim=True)).squeeze()
    return GatedABMIL()


def run_fold(tr_idx, va_idx, cfg):
    n_pca      = cfg['N_PCA']
    hidden_dim = cfg['HIDDEN_DIM']
    att_dim    = cfg['ATT_DIM']
    max_cells  = cfg['MAX_CELLS']
    dropout    = cfg['DROPOUT']
    lr         = cfg['LR']

    pca_feat = pca_full[:, :n_pca]
    pids_tr  = all_pids[tr_idx].tolist()
    pids_va  = all_pids[va_idx].tolist()

    model = build_model(n_pca, hidden_dim, att_dim, dropout)
    opt   = Adam(model.parameters(), lr=lr, weight_decay=WD)

    best_auc, best_state, pat = 0.0, None, 0

    for epoch in range(EPOCHS):
        # train
        model.train()
        order = pids_tr.copy()
        np.random.shuffle(order)
        for pid in order:
            idx = patient_meta[pid]['idx']
            if len(idx) > max_cells:
                idx = np.random.choice(idx, max_cells, replace=False)
            x   = torch.tensor(pca_feat[idx])
            lbl = torch.tensor(patient_meta[pid]['label'], dtype=torch.float32)
            opt.zero_grad()
            F.binary_cross_entropy_with_logits(model(x), lbl).backward()
            opt.step()

        # validate
        model.eval()
        with torch.no_grad():
            probs_va, lbl_va = [], []
            for pid in pids_va:
                idx = patient_meta[pid]['idx']
                probs_va.append(
                    torch.sigmoid(model(torch.tensor(pca_feat[idx]))).item())
                lbl_va.append(patient_meta[pid]['label'])
        probs_va = np.array(probs_va)
        lbl_va   = np.array(lbl_va)

        if len(np.unique(lbl_va)) > 1:
            auc = roc_auc_score(lbl_va, probs_va)
        else:
            auc = 0.5

        if auc > best_auc:
            best_auc   = auc
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            pat = 0
        else:
            pat += 1
        if pat >= PATIENCE:
            break

    # final val with best state
    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        oof_p, oof_l = [], []
        for pid in pids_va:
            idx = patient_meta[pid]['idx']
            oof_p.append(
                torch.sigmoid(model(torch.tensor(pca_feat[idx]))).item())
            oof_l.append(patient_meta[pid]['label'])
    oof_p = np.array(oof_p)
    oof_l = np.array(oof_l)
    auc_final = roc_auc_score(oof_l, oof_p) if len(np.unique(oof_l)) > 1 else 0.5
    return auc_final, oof_p, oof_l


def sweep_param(param_name, values):
    results = []   # list of (value, mean_auc, std_auc, fold_aucs)
    for val in values:
        cfg = dict(BASELINE)
        cfg[param_name] = val
        fold_aucs = []
        oof_p_all = np.zeros(len(all_pids))
        oof_l_all = np.zeros(len(all_pids), dtype=int)
        for tr_idx, va_idx in folds:
            auc, oof_p, oof_l = run_fold(tr_idx, va_idx, cfg)
            fold_aucs.append(auc)
            oof_p_all[va_idx] = oof_p
            oof_l_all[va_idx] = oof_l
        oof_auc = roc_auc_score(oof_l_all, oof_p_all)
        results.append((val, oof_auc, np.std(fold_aucs), fold_aucs))
        print(f'  {param_name}={val}: OOF AUC={oof_auc:.4f}  '
              f'fold={np.mean(fold_aucs):.4f}±{np.std(fold_aucs):.4f}')
    return results


# ══════════════════════════════════════════════════════════════════════════
#  Run sweeps
# ══════════════════════════════════════════════════════════════════════════
all_results = {}
for pname, vals in SWEEPS.items():
    print(f'\n=== Sweeping {pname} ===')
    all_results[pname] = sweep_param(pname, vals)

# ── Save numeric results ────────────────────────────────────────────────────
rows = []
for pname, res in all_results.items():
    for val, oof_auc, fold_std, fold_aucs in res:
        rows.append(dict(
            param=pname, value=val,
            oof_auc=oof_auc, fold_std=fold_std,
            fold_aucs=';'.join(f'{a:.4f}' for a in fold_aucs)
        ))
pd.DataFrame(rows).to_csv('phase3_results/mil/supp_robustness.csv', index=False)
print('\nSaved: phase3_results/mil/supp_robustness.csv')

# ══════════════════════════════════════════════════════════════════════════
#  Plot
# ══════════════════════════════════════════════════════════════════════════
print('\nPlotting Supp Fig 3 ...')

COLOR_LINE    = '#3C5488'
COLOR_BAND    = '#8491B4'
COLOR_BASE    = '#E64B35'   # baseline marker
COLOR_GRID    = '#E8E8E8'

fig, axes = plt.subplots(2, 3, figsize=(11, 7))
axes = axes.flatten()

for ax, (pname, res) in zip(axes, all_results.items()):
    vals      = [r[0] for r in res]
    oof_aucs  = [r[1] for r in res]
    fold_stds = [r[2] for r in res]
    fold_all  = [r[3] for r in res]  # list of lists

    x = np.arange(len(vals))

    # fold range band
    fold_min = [np.min(fa) for fa in fold_all]
    fold_max = [np.max(fa) for fa in fold_all]
    ax.fill_between(x, fold_min, fold_max, alpha=0.15, color=COLOR_LINE,
                    label='Fold range')

    # OOF AUC line
    ax.plot(x, oof_aucs, 'o-', color=COLOR_LINE, lw=1.8, ms=5, zorder=3,
            label='OOF AUC')

    # error bars (fold std)
    ax.errorbar(x, oof_aucs, yerr=fold_stds, fmt='none',
                ecolor=COLOR_BAND, elinewidth=1.2, capsize=3, zorder=4)

    # baseline marker
    base_val = BASELINE[pname]
    if base_val in vals:
        bi = vals.index(base_val)
        ax.axvline(bi, ls='--', color=COLOR_BASE, lw=1.2, alpha=0.7)
        ax.scatter([bi], [oof_aucs[bi]], color=COLOR_BASE, s=55, zorder=5,
                   label=f'Baseline ({base_val})')

    # x-axis ticks
    if pname == 'LR':
        tick_labels = [f'{v:.0e}' for v in vals]
    elif pname == 'DROPOUT':
        tick_labels = [str(v) for v in vals]
    else:
        tick_labels = [str(int(v)) for v in vals]
    ax.set_xticks(x)
    ax.set_xticklabels(tick_labels, fontsize=8.5)

    # y-axis
    ymin_plot = max(0.5, min(fold_min) - 0.05)
    ymax_plot = min(1.0, max(fold_max) + 0.05)
    ax.set_ylim(ymin_plot, ymax_plot)
    ax.set_yticks(np.arange(np.ceil(ymin_plot * 10) / 10,
                             ymax_plot + 0.01, 0.05))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{v:.2f}'))

    ax.set_xlabel(PARAM_LABELS[pname], fontsize=9.5)
    ax.set_ylabel('AUC', fontsize=9.5)
    ax.set_title(pname, fontsize=10.5, fontweight='bold')
    ax.grid(axis='y', color=COLOR_GRID, lw=0.8, zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=7.5, loc='lower right', framealpha=0.85,
              edgecolor='lightgray')

fig.suptitle(
    'Supplementary Figure 3 — ABMIL Hyperparameter Robustness\n'
    'GSE178341 Immune Microenvironment · 5-Fold CV OOF AUC',
    fontsize=11, fontweight='bold', y=1.01
)
plt.tight_layout()

fig.savefig(f'{FIG_DIR}/suppFig3_robustness.png', dpi=200, bbox_inches='tight')
fig.savefig(f'{FIG_DIR}/suppFig3_robustness.pdf', bbox_inches='tight')
plt.close(fig)
print(f'Saved: {FIG_DIR}/suppFig3_robustness.png/.pdf')
print('Done.')
