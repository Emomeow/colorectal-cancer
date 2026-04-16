#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 2-C — TF activity along myeloid pseudotime trajectory
Goal: Identify transcription factors that rise along the
      inflammatory → suppressive myeloid trajectory
      → propose candidate regulators / drug targets
Method: Score TF regulon gene sets in each cell; compute Spearman
        correlation with Palantir pseudotime; identify TFs enriched
        in suppressive (CRC) fate
Depends on: phase2_results/myeloid_v3/ (run phase2_myeloid_v3.py first)
Output: phase2_results/tf_v2/  +  figures/phase2/tf_v2/
"""
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import os, sys, io, warnings
warnings.filterwarnings('ignore')

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42
plt.rcParams['font.family']  = 'Arial'

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

RES = 'phase2_results/tf_v2'
FIG = 'figures/phase2/tf_v2'
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

MYE_RES = 'phase2_results/myeloid_v3'

print('=' * 65)
print('Phase 2-C  TF Trajectory Analysis')
print('=' * 65)

# ─────────────────────────────────────────────────────────────────
# TF regulon gene sets
# Selected from DoRothEA / CancerSEA / literature for myeloid context
# ─────────────────────────────────────────────────────────────────
TF_REGULONS = {
    # NF-κB / inflammatory activation
    'RELA':   ['IL6', 'CXCL8', 'CCL2', 'TNF', 'IL1B', 'PTGS2',
               'ICAM1', 'VCAM1', 'MMP9', 'BIRC2', 'BCL2A1'],
    'NFKB1':  ['IL6', 'CXCL8', 'IL1A', 'TNF', 'PTGS2', 'ICAM1',
               'MMP9', 'BIRC3', 'CCL5', 'CXCL1', 'CXCL2'],
    # Hypoxia / metabolic shift
    'HIF1A':  ['VEGFA', 'LDHA', 'SLC2A1', 'HMOX1', 'BNIP3',
               'CA9', 'ENO1', 'ALDOA', 'GAPDH', 'PGK1', 'ADM'],
    # M2/anti-inflammatory polarization
    'STAT3':  ['IL10', 'SOCS3', 'BCL2', 'MYC', 'VEGFA', 'HIF1A',
               'MMP2', 'TWIST1', 'SNAI1', 'CD44'],
    'PPARG':  ['CD36', 'FABP4', 'ADIPOQ', 'LPL', 'APOC1',
               'MRC1', 'CD163', 'LGMN', 'APOE', 'TREM2'],
    'KLF4':   ['MRC1', 'CD163', 'ARG1', 'CCL22', 'CLEC7A',
               'MGL2', 'RETNLA', 'IL4R', 'TGFB1'],
    'CEBPB':  ['IL10', 'IL6', 'CCL2', 'CXCL8', 'SOCS3',
               'MRC1', 'CD163', 'SPP1', 'FN1', 'LGALS3'],
    # IRFs
    'IRF4':   ['IL10', 'TGFB1', 'CCL17', 'CCL22', 'CD206',
               'CLEC7A', 'CX3CR1', 'MAFB'],
    'IRF1':   ['HLA-DRA', 'HLA-DRB1', 'IFIT1', 'IFIT2', 'ISG15',
               'STAT1', 'OAS1', 'RSAD2', 'MX1', 'CXCL10'],
    # Exhaustion / tolerogenesis
    'NFATC1': ['CD274', 'PDCD1LG2', 'CD80', 'CD86', 'IL10',
               'TGFB1', 'IDO1', 'LGALS9', 'PVR'],
    # AP-1 complex
    'FOS':    ['JUN', 'CCL2', 'CXCL8', 'IL6', 'MMP9',
               'VEGFA', 'FN1', 'FOSL2', 'JUNB'],
    'JUN':    ['FOSL1', 'FOSL2', 'JUNB', 'FOS', 'IL6', 'CCL2',
               'MMP9', 'VEGFC', 'SPP1', 'FN1'],
    # SPP1+ TAM specific
    'MAFB':   ['MRC1', 'APOE', 'C1QA', 'C1QB', 'C1QC', 'SPP1',
               'LGMN', 'FN1', 'CD68', 'TREM2'],
    # c-Myc
    'MYC':    ['CDK4', 'CCND1', 'PCNA', 'MCM2', 'TOP2A',
               'MKI67', 'LDHA', 'FASN', 'PKM'],
    # Runx
    'RUNX1':  ['CEBPA', 'SPI1', 'MYB', 'ITGAM', 'CSF1R',
               'FLT3', 'GATA2', 'CBFB'],
    # SPI1 / PU.1 (myeloid master TF)
    'SPI1':   ['CSF1R', 'ITGAM', 'CD68', 'MRC1', 'FCGR3A',
               'CYBB', 'NCF1', 'CEBPA', 'IRF8'],
}

print(f'  TFs to score: {len(TF_REGULONS)}')

# ─────────────────────────────────────────────────────────────────
# 1. Load myeloid cells
# ─────────────────────────────────────────────────────────────────
print('\n[1/5] Loading myeloid data ...')
adata_all = sc.read_h5ad('processed_data/merged_annotated_fine.h5ad')

mask = (adata_all.obs['batch'].isin(['GSE125527', 'GSE132465'])) & \
       (adata_all.obs['cell_type_grouped'] == 'Myeloids')
adata = adata_all[mask].copy()
del adata_all
print(f'  Cells: {adata.n_obs:,}')

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Load pseudotime
pt_csv = f'{MYE_RES}/pseudotime_values.csv'
if not os.path.exists(pt_csv):
    sys.exit(f'ERROR: {pt_csv} not found. Run phase2_myeloid_v3.py first.')
pt_vals = pd.read_csv(pt_csv, index_col=0).squeeze().reindex(adata.obs_names)
adata.obs['pseudotime'] = pt_vals.values

# Load branch probs
bp_csv = f'{MYE_RES}/branch_probs.csv'
bp_df  = pd.read_csv(bp_csv, index_col=0).reindex(adata.obs_names)
adata.obs['is_predriver'] = pd.read_csv(
    f'{MYE_RES}/cell_metadata.csv', index_col=0
).reindex(adata.obs_names)['supp_score'].values > \
    pd.read_csv(f'{MYE_RES}/cell_metadata.csv', index_col=0)['supp_score'].median()

print(f'  Pseudotime: {pt_vals.min():.3f} – {pt_vals.max():.3f}')
print(f'  NaN pseudotime: {pt_vals.isna().sum()}')

# ─────────────────────────────────────────────────────────────────
# 2. Score TF regulons
# ─────────────────────────────────────────────────────────────────
print('\n[2/5] Scoring TF regulons ...')
avail_genes = set(adata.var_names)
tf_stats = {}

for tf, gene_set in TF_REGULONS.items():
    use_genes = [g for g in gene_set if g in avail_genes]
    if len(use_genes) >= 3:
        sc.tl.score_genes(adata, use_genes, score_name=f'tf_{tf}', use_raw=False)
        tf_stats[tf] = {'n_genes': len(use_genes), 'total': len(gene_set)}
        print(f'  {tf}: {len(use_genes)}/{len(gene_set)} genes')
    else:
        print(f'  {tf}: only {len(use_genes)} genes — skipped')

tf_scored = list(tf_stats.keys())
print(f'\n  Scored {len(tf_scored)} TFs')

# ─────────────────────────────────────────────────────────────────
# 3. Spearman correlation with pseudotime + fate DEG
# ─────────────────────────────────────────────────────────────────
print('\n[3/5] Correlating TF scores with pseudotime ...')
valid = ~pt_vals.isna().values
pt_valid = pt_vals.values[valid]

corr_rows = []
for tf in tf_scored:
    col = f'tf_{tf}'
    scores = adata.obs[col].values[valid]
    rho, pval = stats.spearmanr(scores, pt_valid)
    corr_rows.append({
        'TF': tf,
        'rho': rho,
        'pval': pval,
        'abs_rho': abs(rho),
        'direction': 'up' if rho > 0 else 'down',
    })

corr_df = pd.DataFrame(corr_rows).sort_values('abs_rho', ascending=False)
corr_df.to_csv(f'{RES}/tf_pseudotime_correlation.csv', index=False)
print(corr_df.to_string(index=False))

# Also compute mean TF activity per subtype
subtype_tf = {}
for st in adata.obs['cell_subtype_fine'].unique():
    m = adata.obs['cell_subtype_fine'] == st
    row = {tf: adata.obs.loc[m, f'tf_{tf}'].mean() for tf in tf_scored}
    row['subtype'] = st
    subtype_tf[st] = row
subtype_tf_df = pd.DataFrame(subtype_tf).T
subtype_tf_df.to_csv(f'{RES}/tf_activity_by_subtype.csv')
print(f'\n  TF activity by subtype:')
print(subtype_tf_df[corr_df['TF'].head(8).tolist()].to_string())

# ─────────────────────────────────────────────────────────────────
# 4. TF trend along pseudotime quantiles
# ─────────────────────────────────────────────────────────────────
print('\n[4/5] TF trends along pseudotime ...')
N_BINS = 20
pt_bins = pd.qcut(pt_vals.values, q=N_BINS, labels=False, duplicates='drop')
adata.obs['pt_bin'] = pt_bins

trend_rows = []
for b in range(N_BINS):
    m = adata.obs['pt_bin'] == b
    if m.sum() < 5:
        continue
    row = {'pt_bin': b,
           'pt_mean': pt_vals.values[m].mean(),
           'n_cells': m.sum()}
    for tf in tf_scored:
        row[tf] = adata.obs.loc[m, f'tf_{tf}'].mean()
    trend_rows.append(row)

trend_df = pd.DataFrame(trend_rows).set_index('pt_bin')
trend_df.to_csv(f'{RES}/tf_trends_along_pseudotime.csv')

# ─────────────────────────────────────────────────────────────────
# 5. Plots
# ─────────────────────────────────────────────────────────────────
print('\n[5/5] Plotting ...')

# ── Fig A: TF correlation ranking ─────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
colors = ['#E64B35' if r > 0 else '#4DBBD5' for r in corr_df['rho']]
ax.barh(range(len(corr_df)), corr_df['rho'].values, color=colors, alpha=0.85)
ax.set_yticks(range(len(corr_df)))
ax.set_yticklabels(corr_df['TF'].values, fontsize=9, style='italic')
ax.axvline(0, color='black', linewidth=0.8)
ax.set_xlabel("Spearman ρ (TF score vs pseudotime)", fontsize=10)
ax.set_title("TF activity correlation with myeloid pseudotime", fontsize=11, fontweight='bold')
red_patch = mpatches.Patch(color='#E64B35', label='Rising (suppressive)')
blue_patch = mpatches.Patch(color='#4DBBD5', label='Falling (inflammatory)')
ax.legend(handles=[red_patch, blue_patch], fontsize=9)
fig.tight_layout()
for ext in ('png', 'pdf'):
    path = f'{FIG}/tf_pseudotime_correlation.{ext}'
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f'  Saved: {path}')
plt.close()

# ── Fig B: TF heatmap (subtypes × TFs) ────────────────────────
SUBTYPE_ORDER = ['Classical_Mono', 'Inflam_Mac', 'TAM', 'cDC2', 'pDC']
plot_tf  = [t for t in corr_df['TF'].head(12).tolist() if t in subtype_tf_df.columns]
plot_sub = [s for s in SUBTYPE_ORDER if s in subtype_tf_df.index]

if plot_sub and plot_tf:
    mat = subtype_tf_df.loc[plot_sub, plot_tf].values.astype(float)
    # Z-score per TF
    mat_z = (mat - mat.mean(axis=0)) / (mat.std(axis=0) + 1e-8)
    mat_z = np.clip(mat_z, -2.5, 2.5)

    fig, ax = plt.subplots(figsize=(len(plot_tf) * 0.85 + 1.5, len(plot_sub) * 0.7 + 1.5))
    im = ax.imshow(mat_z, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(len(plot_tf)))
    ax.set_xticklabels(plot_tf, rotation=45, ha='right', fontsize=9, style='italic')
    ax.set_yticks(range(len(plot_sub)))
    ax.set_yticklabels(plot_sub, fontsize=9)
    plt.colorbar(im, ax=ax, shrink=0.8, label='Z-score')
    ax.set_title('TF regulon activity by myeloid subtype', fontsize=11, fontweight='bold')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        path = f'{FIG}/tf_heatmap_subtype.{ext}'
        fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
        print(f'  Saved: {path}')
    plt.close()

# ── Fig C: TF trend curves along pseudotime ───────────────────
# Top rising TFs (positive correlation)
rising_tfs = corr_df[corr_df['rho'] > 0]['TF'].head(6).tolist()
falling_tfs = corr_df[corr_df['rho'] < 0]['TF'].head(3).tolist()
plot_tfs = rising_tfs + falling_tfs
plot_tfs = [t for t in plot_tfs if t in trend_df.columns]

if plot_tfs and len(trend_df) > 3:
    n = len(plot_tfs)
    nc = 3
    nr = (n + nc - 1) // nc
    fig, axes = plt.subplots(nr, nc, figsize=(5 * nc, 3.5 * nr))
    axes = np.array(axes).ravel()
    fig.suptitle('TF regulon activity along myeloid pseudotime',
                 fontsize=11, fontweight='bold', fontfamily='Arial')

    COLORS_TF = ['#E64B35', '#F39B7F', '#DC0000', '#B22222', '#FF6347', '#FF4500',
                 '#4DBBD5', '#3C5488', '#00A087']

    for i, (ax, tf) in enumerate(zip(axes, plot_tfs)):
        rho = corr_df.loc[corr_df['TF'] == tf, 'rho'].values[0]
        col = COLORS_TF[i % len(COLORS_TF)]
        ax.plot(trend_df['pt_mean'], trend_df[tf], color=col, linewidth=2)
        ax.fill_between(trend_df['pt_mean'], trend_df[tf], alpha=0.2, color=col)
        ax.set_xlabel('Pseudotime', fontsize=8)
        ax.set_ylabel('Mean TF score', fontsize=8)
        ax.set_title(f'{tf}  (ρ={rho:+.2f})', fontsize=10, fontweight='bold',
                     style='italic', color=col)
        ax.tick_params(labelsize=7)

    for j in range(i + 1, len(axes)):
        axes[j].axis('off')

    fig.tight_layout()
    for ext in ('png', 'pdf'):
        path = f'{FIG}/tf_trends.{ext}'
        fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
        print(f'  Saved: {path}')
    plt.close()

# ── Fig D: Top TFs on UMAP ────────────────────────────────────
# Reload UMAP from myeloid_v3 metadata
# Recompute neighbors + UMAP for visualization
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=30, random_state=0)
sc.tl.umap(adata, min_dist=0.3, random_state=0)

umap = adata.obsm['X_umap']
shuf = np.random.default_rng(0).permutation(adata.n_obs)
top_tfs = corr_df['TF'].head(6).tolist()
top_tfs = [t for t in top_tfs if f'tf_{t}' in adata.obs.columns]

if top_tfs:
    nc = 3
    nr = (len(top_tfs) + nc - 1) // nc
    fig, axes = plt.subplots(nr, nc, figsize=(5 * nc, 4 * nr))
    axes = np.array(axes).ravel()
    fig.suptitle('Top TF regulon scores on myeloid UMAP',
                 fontsize=11, fontweight='bold', fontfamily='Arial')

    for i, (ax, tf) in enumerate(zip(axes, top_tfs)):
        scores = adata.obs[f'tf_{tf}'].values
        sc_ = ax.scatter(umap[shuf, 0], umap[shuf, 1],
                         c=scores[shuf], s=2, alpha=0.7,
                         cmap='RdYlBu_r', linewidths=0)
        plt.colorbar(sc_, ax=ax, shrink=0.8, pad=0.01)
        rho = corr_df.loc[corr_df['TF'] == tf, 'rho'].values[0]
        ax.set_title(f'{tf} (ρ={rho:+.2f})', fontsize=10, fontweight='bold', style='italic')
        ax.axis('off')

    for j in range(i + 1, len(axes)):
        axes[j].axis('off')

    fig.tight_layout()
    for ext in ('png', 'pdf'):
        path = f'{FIG}/tf_umap.{ext}'
        fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
        print(f'  Saved: {path}')
    plt.close()

print('\n' + '=' * 65)
print('Summary — TF Trajectory (Phase 2-C)')
print('=' * 65)
print(f'  TFs scored: {len(tf_scored)}')
print(f'\n  Top rising TFs (suppressive branch):')
print(corr_df[corr_df['rho'] > 0].head(8)[['TF', 'rho', 'pval']].to_string(index=False))
print(f'\n  Top falling TFs (inflammatory):')
print(corr_df[corr_df['rho'] < 0].head(5)[['TF', 'rho', 'pval']].to_string(index=False))
print(f'\n  Results in : {RES}/')
print(f'  Figures in : {FIG}/')
