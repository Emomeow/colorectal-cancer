#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 2-A — Myeloid Pseudotime v3
Data   : GSE125527 (UC, inflammatory) + GSE132465 (CRC, tumor)  [ONLY these two]
Goal   : Trajectory from inflammatory monocyte → suppressive TAM/Pre-driver
Root   : Classical_Mono from UC with lowest suppressive gene score
Output : phase2_results/myeloid_v3/  +  figures/phase2/myeloid_v3/
"""
import scanpy as sc
import numpy as np
import pandas as pd
import palantir
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import os, sys, io, warnings
from datetime import datetime
warnings.filterwarnings('ignore')

# ── font / PDF settings ───────────────────────────────────────────
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42
plt.rcParams['font.family']  = 'Arial'

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')


def find_branch_point(branch_probs, pseudotime, n_bins=50):
    """
    Branch point = pseudotime at which fate-probability variance is maximal.
    This is the data-driven transition point where fate commitment begins.
    """
    pt = pseudotime.values if hasattr(pseudotime, 'values') else pseudotime
    valid = ~np.isnan(pt)
    pt_v, bp_v = pt[valid], branch_probs.values[valid]
    bins    = np.linspace(pt_v.min(), pt_v.max(), n_bins + 1)
    bin_idx = np.clip(np.digitize(pt_v, bins) - 1, 0, n_bins - 1)
    variances = [bp_v[bin_idx == b].var(axis=0).mean()
                 if (bin_idx == b).sum() > 2 else 0
                 for b in range(n_bins)]
    peak = np.argmax(variances)
    return ((bins[:-1] + bins[1:]) / 2)[peak]

RES = 'phase2_results/myeloid_v3'
FIG = 'figures/phase2/myeloid_v3'
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

print('=' * 65)
print(f'Phase 2-A  Myeloid Pseudotime v3  |  {datetime.now().strftime("%H:%M:%S")}')
print('=' * 65)

# ─────────────────────────────────────────────────────────────────
# 1. Load & subset
# ─────────────────────────────────────────────────────────────────
print('\n[1/9] Loading data ...')
adata_all = sc.read_h5ad('processed_data/merged_annotated_fine.h5ad', backed='r')

mask = (adata_all.obs['batch'].isin(['GSE125527', 'GSE132465'])) & \
       (adata_all.obs['cell_type_grouped'] == 'Myeloids')
adata = adata_all[mask].to_memory()
adata_all.file.close()
del adata_all

print(f'  Subset: {adata.n_obs:,} cells')
print(pd.crosstab(adata.obs['batch'], adata.obs['cell_subtype_fine']).to_string())

# Subtype color palette
SUBTYPE_COLORS = {
    'Classical_Mono': '#4DBBD5',
    'Inflam_Mac'    : '#E64B35',
    'TAM'           : '#F39B7F',
    'cDC2'          : '#00A087',
    'pDC'           : '#8491B4',
    'Unknown_My'    : '#B8B8B8',
}
BATCH_COLORS = {'GSE125527': '#3C5488', 'GSE132465': '#DC0000'}

# ─────────────────────────────────────────────────────────────────
# 2. Normalize expression (for gene scoring & DEG)
# ─────────────────────────────────────────────────────────────────
print('\n[2/9] Normalizing expression ...')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ─────────────────────────────────────────────────────────────────
# 3. KNN + UMAP on X_scvi of this subset
# ─────────────────────────────────────────────────────────────────
print('\n[3/9] KNN + UMAP on X_scvi subset ...')
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=30, random_state=0)
sc.tl.umap(adata, min_dist=0.3, random_state=0)
print(f'  UMAP: {adata.obsm["X_umap"].shape}')

# ─────────────────────────────────────────────────────────────────
# 4. Biologically-motivated root selection
#    Root = Classical_Mono from UC with lowest suppressive score
# ─────────────────────────────────────────────────────────────────
print('\n[4/9] Selecting root cell ...')

# Suppressive gene module (inflammatory → suppressive trajectory markers)
SUPP_GENES = ['SPP1', 'APOE', 'C1QC', 'IL10', 'TGFB1', 'LGALS3', 'MRC1', 'VEGFA', 'CD274']
avail_supp = [g for g in SUPP_GENES if g in adata.var_names]
print(f'  Suppressive genes available: {avail_supp}')

# Score suppressive module
sc.tl.score_genes(adata, avail_supp, score_name='supp_score', use_raw=False)

# Root candidates: Classical_Mono from UC (GSE125527)
cand_mask = (adata.obs['batch'] == 'GSE125527') & \
            (adata.obs['cell_subtype_fine'] == 'Classical_Mono')
cand_idx  = np.where(cand_mask.values)[0]
supp_vals = adata.obs['supp_score'].values[cand_idx]
root_local = cand_idx[np.argmin(supp_vals)]
root_cell  = adata.obs_names[root_local]

print(f'  Root: {root_cell}')
print(f'    subtype={adata.obs.loc[root_cell, "cell_subtype_fine"]}  '
      f'batch={adata.obs.loc[root_cell, "batch"]}  '
      f'supp_score={adata.obs.loc[root_cell, "supp_score"]:.3f}')

# ─────────────────────────────────────────────────────────────────
# 5. Diffusion map + Palantir
# ─────────────────────────────────────────────────────────────────
print('\n[5/9] Diffusion map ...')
dm_res = palantir.utils.run_diffusion_maps(
    pd.DataFrame(adata.obsm['X_scvi'], index=adata.obs_names),
    n_components=10, knn=30,
)
adata.obsm['X_diffmap'] = dm_res['EigenVectors'].values
ms_data = palantir.utils.determine_multiscale_space(dm_res)

print('\n[5b/9] Palantir pseudotime ...')
pr_res = palantir.core.run_palantir(
    ms_data, early_cell=root_cell,
    num_waypoints=1200, use_early_cell_as_start=True, n_jobs=4,
)
adata.obs['palantir_pseudotime'] = pr_res.pseudotime
adata.obs['entropy']             = pr_res.entropy
terminal_ids = list(pr_res.branch_probs.columns)
n_term       = len(terminal_ids)

print(f'  Pseudotime: {pr_res.pseudotime.min():.3f} – {pr_res.pseudotime.max():.3f}')
print(f'  Terminals ({n_term}):')
for t in terminal_ids:
    print(f'    {t}  subtype={adata.obs.loc[t, "cell_subtype_fine"]}  '
          f'batch={adata.obs.loc[t, "batch"]}')

# Save core results
pr_res.pseudotime.to_csv(f'{RES}/pseudotime_values.csv')
pr_res.branch_probs.to_csv(f'{RES}/branch_probs.csv')
adata.obs[['batch', 'cell_subtype_fine', 'palantir_pseudotime',
           'entropy', 'supp_score']].to_csv(f'{RES}/cell_metadata.csv')

# ─────────────────────────────────────────────────────────────────
# 6. Define Pre-driver candidates
#    Criteria (from Phase2 mission):
#      a) pseudotime within branch-point window (data-driven)
#         window = [branch_pt - 0.30 * span, branch_pt]
#         where branch_pt = pseudotime of max fate-prob variance
#      b) supp_score above median of Classical_Mono
#      c) branch probability toward TAM/suppressive terminal > 0.5
#      d) present in BOTH GSE125527 and GSE132465
# ─────────────────────────────────────────────────────────────────
print('\n[6/9] Defining Pre-driver candidates ...')

pt = adata.obs['palantir_pseudotime'].values

# Branch point-based window (replaces fixed [40th, 80th] percentile)
branch_pt   = find_branch_point(pr_res.branch_probs, pr_res.pseudotime)
pt_span     = pt.max() - pt.min()
pt_lo       = branch_pt - 0.30 * pt_span   # 30% span before branch point
pt_hi       = branch_pt                    # upper boundary = branch point
print(f'  Branch point pseudotime: {branch_pt:.3f}')
print(f'  Pre-driver window: [{pt_lo:.3f}, {pt_hi:.3f}]'
      f'  (was fixed [40th,80th]={np.percentile(pt,40):.3f},{np.percentile(pt,80):.3f})')

# Baseline suppressive score threshold (median of Classical_Mono)
mono_supp_median = adata.obs.loc[
    adata.obs['cell_subtype_fine'] == 'Classical_Mono', 'supp_score'
].median()

# Identify "suppressive terminal" = terminal in GSE132465 (CRC)
bp = pr_res.branch_probs.reindex(adata.obs_names)
crc_terms = [t for t in terminal_ids
             if adata.obs.loc[t, 'batch'] == 'GSE132465']
uc_terms  = [t for t in terminal_ids
             if adata.obs.loc[t, 'batch'] == 'GSE125527']
print(f'  CRC terminals: {crc_terms}')
print(f'  UC  terminals: {uc_terms}')

p_supp = bp[crc_terms].sum(axis=1).values if crc_terms else np.zeros(adata.n_obs)
p_infl = bp[uc_terms].sum(axis=1).values  if uc_terms  else np.zeros(adata.n_obs)

# Pre-driver mask
crit_a = (pt >= pt_lo) & (pt <= pt_hi)
crit_b = adata.obs['supp_score'].values > mono_supp_median
crit_c = p_supp > 0.5
predriver_mask = crit_a & crit_b & crit_c
adata.obs['is_predriver'] = predriver_mask

n_pd = predriver_mask.sum()
print(f'  Pre-driver candidates: {n_pd:,} cells  '
      f'({n_pd / adata.n_obs * 100:.1f}%)')
print('  Subtype breakdown:')
print(adata.obs.loc[predriver_mask, 'cell_subtype_fine'].value_counts().to_string())
print('  Batch breakdown:')
print(adata.obs.loc[predriver_mask, 'batch'].value_counts().to_string())

# Save
predriver_df = adata.obs[predriver_mask][['batch', 'cell_subtype_fine',
                                          'palantir_pseudotime', 'supp_score']].copy()
predriver_df['p_supp'] = p_supp[predriver_mask]
predriver_df.to_csv(f'{RES}/predriver_candidates.csv')
print(f'  Saved: {RES}/predriver_candidates.csv')

# ─────────────────────────────────────────────────────────────────
# 7. Branch genes (Spearman correlation + fate DEG)
# ─────────────────────────────────────────────────────────────────
print('\n[7/9] Branch gene analysis ...')

# 7a. Spearman correlation with pseudotime (all HVGs)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='batch')
hvg_genes = adata.var_names[adata.var['highly_variable']].tolist()
X_hvg = adata[:, hvg_genes].X
if hasattr(X_hvg, 'toarray'):
    X_hvg = X_hvg.toarray()

pt_vals = adata.obs['palantir_pseudotime'].values
rho_list = []
for i in range(X_hvg.shape[1]):
    rho, pval = stats.spearmanr(X_hvg[:, i], pt_vals)
    rho_list.append((hvg_genes[i], rho, pval))

pt_corr = pd.DataFrame(rho_list, columns=['gene', 'rho', 'pval'])
pt_corr['abs_rho'] = pt_corr['rho'].abs()
pt_corr = pt_corr.sort_values('abs_rho', ascending=False)
pt_corr.to_csv(f'{RES}/pseudotime_gene_correlation.csv', index=False)
print(f'  Top genes correlated with pseudotime:')
print(pt_corr.head(20).to_string(index=False))

# 7b. Fate DEG (CRC-fate vs UC-fate cells)
deg_crc = deg_uc = None
if crc_terms and uc_terms:
    thr_crc = np.percentile(p_supp, 70)
    thr_uc  = np.percentile(p_infl, 70)
    keep    = (p_supp > thr_crc) | (p_infl > thr_uc)
    adata_bp = adata[keep].copy()
    adata_bp.obs['fate_label'] = np.where(
        p_supp[keep] > p_infl[keep], 'CRC_fate', 'UC_fate'
    )
    n_crc = (adata_bp.obs['fate_label'] == 'CRC_fate').sum()
    n_uc  = (adata_bp.obs['fate_label'] == 'UC_fate').sum()
    print(f'\n  Fate DEG: CRC_fate={n_crc:,}  UC_fate={n_uc:,}')

    for k in ['neighbors', 'diffmap_evals']:
        adata_bp.uns.pop(k, None)
    sc.tl.rank_genes_groups(
        adata_bp, groupby='fate_label',
        groups=['CRC_fate', 'UC_fate'], reference='rest',
        method='wilcoxon', key_added='branch_deg',
    )
    deg_crc = sc.get.rank_genes_groups_df(adata_bp, group='CRC_fate', key='branch_deg')
    deg_uc  = sc.get.rank_genes_groups_df(adata_bp, group='UC_fate',  key='branch_deg')
    deg_crc = deg_crc[deg_crc['logfoldchanges'] > 0].head(50)
    deg_uc  = deg_uc[deg_uc['logfoldchanges'] > 0].head(50)
    deg_crc.to_csv(f'{RES}/branch_DEG_CRC_fate.csv', index=False)
    deg_uc.to_csv(f'{RES}/branch_DEG_UC_fate.csv',   index=False)
    print(f'  Top CRC_fate genes: {", ".join(deg_crc["names"].head(10))}')
    print(f'  Top UC_fate  genes: {", ".join(deg_uc["names"].head(10))}')
elif n_term >= 2:
    # Fallback: compare any two terminals
    thr1 = np.percentile(bp.iloc[:, 0].values, 70)
    thr2 = np.percentile(bp.iloc[:, 1].values, 70)
    keep = (bp.iloc[:, 0].values > thr1) | (bp.iloc[:, 1].values > thr2)
    adata_bp = adata[keep].copy()
    adata_bp.obs['fate_label'] = np.where(
        bp.iloc[:, 0].values[keep] > bp.iloc[:, 1].values[keep], 'fate1', 'fate2'
    )
    for k in ['neighbors', 'diffmap_evals']:
        adata_bp.uns.pop(k, None)
    sc.tl.rank_genes_groups(
        adata_bp, groupby='fate_label',
        groups=['fate1', 'fate2'], reference='rest',
        method='wilcoxon', key_added='branch_deg',
    )
    deg_crc = sc.get.rank_genes_groups_df(adata_bp, group='fate1', key='branch_deg')
    deg_uc  = sc.get.rank_genes_groups_df(adata_bp, group='fate2', key='branch_deg')
    deg_crc = deg_crc[deg_crc['logfoldchanges'] > 0].head(50)
    deg_uc  = deg_uc[deg_uc['logfoldchanges'] > 0].head(50)
    deg_crc.to_csv(f'{RES}/branch_DEG_fate1.csv', index=False)
    deg_uc.to_csv(f'{RES}/branch_DEG_fate2.csv',  index=False)

# ─────────────────────────────────────────────────────────────────
# 8. Gene trends for key suppressive genes
# ─────────────────────────────────────────────────────────────────
print('\n[8/9] Gene trends ...')
KEY_GENES = ['SPP1', 'APOE', 'C1QC', 'IL10', 'TGFB1', 'CD274',
             'MRC1', 'LGALS3', 'VEGFA', 'FN1', 'CCL2', 'CXCL8',
             'HLA-DRA', 'HLA-DRB1', 'FCGR3A', 'S100A8', 'S100A9',
             'HMOX1', 'IDO1', 'PDCD1LG2']

adata_hvg = adata[:, adata.var['highly_variable']].copy()
adata_hvg.obs['palantir_pseudotime'] = adata.obs['palantir_pseudotime']
adata_hvg.obsm['palantir_fate_probabilities'] = bp

try:
    palantir.presults.select_branch_cells(adata_hvg, q=0.01, eps=0.01)
    gene_trends = palantir.presults.compute_gene_trends(
        adata_hvg, lineages=terminal_ids,
        pseudo_time_key='palantir_pseudotime',
    )
    print(f'  Gene trends computed for {len(gene_trends)} lineages')

    # Save gene trend values for key suppressive genes
    for term_id in terminal_ids:
        if term_id in gene_trends:
            trends_df = gene_trends[term_id]['trends']
            batch_lbl = adata.obs.loc[term_id, 'batch']
            avail_key = [g for g in KEY_GENES if g in trends_df.index]
            if avail_key:
                trends_df.loc[avail_key].to_csv(
                    f'{RES}/gene_trends_{batch_lbl}_terminal.csv'
                )
except Exception as e:
    print(f'  Gene trends failed: {e}')
    gene_trends = {}

# ─────────────────────────────────────────────────────────────────
# 9. Plots
# ─────────────────────────────────────────────────────────────────
print('\n[9/9] Plotting ...')
umap  = adata.obsm['X_umap']
pt    = adata.obs['palantir_pseudotime'].values
ent   = adata.obs['entropy'].values
shuf  = np.random.default_rng(0).permutation(adata.n_obs)
S, A  = 2.0, 0.6

# ── Fig A: Overview (4-panel) ──────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Myeloid Pseudotime — GSE125527 (UC) + GSE132465 (CRC)',
             fontsize=12, fontweight='bold', fontfamily='Arial')

# A1: By batch
ax = axes[0, 0]
for bat, col in BATCH_COLORS.items():
    m = adata.obs['batch'].values == bat
    ax.scatter(umap[m, 0], umap[m, 1], c=col, s=S, alpha=A,
               linewidths=0, label=bat)
ax.set_title('Dataset', fontsize=10, fontweight='bold')
ax.axis('off')
ax.legend(markerscale=5, fontsize=9)

# A2: By subtype
ax = axes[0, 1]
for st, col in SUBTYPE_COLORS.items():
    m = adata.obs['cell_subtype_fine'].values == st
    ax.scatter(umap[m, 0], umap[m, 1], c=col, s=S, alpha=A,
               linewidths=0, label=st)
ax.set_title('Cell subtype', fontsize=10, fontweight='bold')
ax.axis('off')
ax.legend(markerscale=5, fontsize=8)

# A3: Pseudotime
ax = axes[1, 0]
sc_pt = ax.scatter(umap[shuf, 0], umap[shuf, 1], c=pt[shuf],
                   s=S, alpha=0.7, cmap='viridis', linewidths=0)
plt.colorbar(sc_pt, ax=ax, shrink=0.8, label='Pseudotime')
root_idx = adata.obs_names.get_loc(root_cell)
ax.scatter(umap[root_idx, 0], umap[root_idx, 1],
           s=250, c='red', marker='*', zorder=10, label='Root')
for t in terminal_ids:
    ti = adata.obs_names.get_loc(t)
    ax.scatter(umap[ti, 0], umap[ti, 1],
               s=150, c='orange', marker='D', zorder=10)
ax.set_title('Pseudotime (* root  ◆ terminals)', fontsize=10, fontweight='bold')
ax.axis('off')
ax.legend(fontsize=9)

# A4: Pre-driver candidates
ax = axes[1, 1]
ax.scatter(umap[~predriver_mask, 0], umap[~predriver_mask, 1],
           c='#CCCCCC', s=S * 0.7, alpha=0.4, linewidths=0)
ax.scatter(umap[predriver_mask, 0], umap[predriver_mask, 1],
           c='#E64B35', s=S * 2, alpha=0.8, linewidths=0, label=f'Pre-driver (n={n_pd:,})')
ax.set_title('Pre-driver candidates', fontsize=10, fontweight='bold')
ax.axis('off')
ax.legend(markerscale=3, fontsize=9)

fig.tight_layout()
for ext in ('png', 'pdf'):
    path = f'{FIG}/pseudotime_overview.{ext}'
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f'  Saved: {path}')
plt.close()

# ── Fig B: Suppressive score + pseudotime distribution ────────
fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
fig.suptitle('Myeloid trajectory — suppressive score & pseudotime',
             fontsize=11, fontweight='bold', fontfamily='Arial')

# B1: Suppressive score on UMAP
ax = axes[0]
sc_ss = ax.scatter(umap[shuf, 0], umap[shuf, 1],
                   c=adata.obs['supp_score'].values[shuf],
                   s=S, alpha=0.7, cmap='YlOrRd', linewidths=0)
plt.colorbar(sc_ss, ax=ax, shrink=0.8, label='Suppressive score')
ax.set_title('Suppressive gene score', fontsize=10, fontweight='bold')
ax.axis('off')

# B2: Pseudotime by subtype (violin)
ax = axes[1]
subtypes_order = ['Classical_Mono', 'Inflam_Mac', 'TAM', 'cDC2', 'pDC']
pt_by_st = [adata.obs.loc[adata.obs['cell_subtype_fine'] == st,
                           'palantir_pseudotime'].values
            for st in subtypes_order]
bp_obj = ax.boxplot(pt_by_st, patch_artist=True, notch=False,
                    labels=subtypes_order, widths=0.6)
for patch, st in zip(bp_obj['boxes'], subtypes_order):
    patch.set_facecolor(SUBTYPE_COLORS.get(st, '#AAAAAA'))
    patch.set_alpha(0.8)
ax.set_ylabel('Pseudotime', fontsize=9)
ax.set_title('Pseudotime by subtype', fontsize=10, fontweight='bold')
ax.tick_params(axis='x', labelsize=8, labelrotation=30)

# B3: Pre-driver vs non-Pre-driver pseudotime distribution
ax = axes[2]
ax.hist(pt[~predriver_mask], bins=40, color='#CCCCCC', alpha=0.7,
        label='Other', density=True)
ax.hist(pt[predriver_mask], bins=40, color='#E64B35', alpha=0.8,
        label=f'Pre-driver (n={n_pd})', density=True)
ax.axvline(pt_lo, color='gray', linestyle='--', linewidth=1)
ax.axvline(pt_hi, color='gray', linestyle='--', linewidth=1)
ax.set_xlabel('Pseudotime', fontsize=9)
ax.set_ylabel('Density', fontsize=9)
ax.set_title('Pre-driver pseudotime window', fontsize=10, fontweight='bold')
ax.legend(fontsize=8)

fig.tight_layout()
for ext in ('png', 'pdf'):
    path = f'{FIG}/pseudotime_distribution.{ext}'
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f'  Saved: {path}')
plt.close()

# ── Fig C: Key gene expression on UMAP ────────────────────────
genes_plot = [g for g in ['SPP1', 'APOE', 'C1QC', 'MRC1',
                           'S100A8', 'HLA-DRA', 'CD274', 'VEGFA']
              if g in adata.var_names]
nc = min(4, len(genes_plot))
nr = (len(genes_plot) + nc - 1) // nc
fig, axes = plt.subplots(nr, nc, figsize=(4 * nc, 3.5 * nr))
axes = np.array(axes).ravel()
fig.suptitle('Key gene expression — Myeloid UMAP',
             fontsize=11, fontweight='bold', fontfamily='Arial')

for i, gene in enumerate(genes_plot):
    ax = axes[i]
    expr = adata[:, gene].X
    if hasattr(expr, 'toarray'):
        expr = expr.toarray().ravel()
    else:
        expr = np.array(expr).ravel()
    sc_ = ax.scatter(umap[shuf, 0], umap[shuf, 1],
                     c=expr[shuf], s=S, alpha=0.7,
                     cmap='RdYlBu_r', linewidths=0)
    plt.colorbar(sc_, ax=ax, shrink=0.8, pad=0.01)
    ax.set_title(gene, fontsize=10, fontweight='bold', style='italic')
    ax.axis('off')

for j in range(len(genes_plot), len(axes)):
    axes[j].axis('off')

fig.tight_layout()
for ext in ('png', 'pdf'):
    path = f'{FIG}/gene_expression_umap.{ext}'
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f'  Saved: {path}')
plt.close()

# ── Fig D: Top branch genes volcano ───────────────────────────
if deg_crc is not None and deg_uc is not None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('Branch DEG: CRC fate vs UC fate (myeloid)',
                 fontsize=11, fontweight='bold', fontfamily='Arial')
    for ax, df, col, label in [
        (axes[0], deg_crc, '#DC0000', 'CRC fate'),
        (axes[1], deg_uc,  '#3C5488', 'UC fate'),
    ]:
        ax.bar(range(min(20, len(df))), df['scores'].head(20), color=col, alpha=0.8)
        ax.set_xticks(range(min(20, len(df))))
        ax.set_xticklabels(df['names'].head(20), rotation=45, ha='right',
                           fontsize=8, style='italic')
        ax.set_ylabel('Wilcoxon score', fontsize=9)
        ax.set_title(f'Top {label} genes', fontsize=10, fontweight='bold')

    fig.tight_layout()
    for ext in ('png', 'pdf'):
        path = f'{FIG}/branch_DEG.{ext}'
        fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
        print(f'  Saved: {path}')
    plt.close()

# ── Fig E: Pseudotime vs supp_score scatter ───────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
fig.suptitle('Pre-driver trajectory summary',
             fontsize=11, fontweight='bold', fontfamily='Arial')

ax = axes[0]
for st, col in SUBTYPE_COLORS.items():
    m = adata.obs['cell_subtype_fine'].values == st
    ax.scatter(pt[m], adata.obs['supp_score'].values[m],
               c=col, s=1.5, alpha=0.5, linewidths=0, label=st)
ax.axvspan(pt_lo, pt_hi, alpha=0.12, color='#E64B35', label='Pre-driver window')
ax.set_xlabel('Pseudotime', fontsize=9)
ax.set_ylabel('Suppressive score', fontsize=9)
ax.set_title('Pseudotime × Suppressive score', fontsize=10, fontweight='bold')
handles = [mpatches.Patch(color=c, label=s)
           for s, c in SUBTYPE_COLORS.items() if s != 'Unknown_My']
handles.append(mpatches.Patch(color='#E64B35', alpha=0.3, label='Pre-driver window'))
ax.legend(handles=handles, fontsize=7, markerscale=3)

ax = axes[1]
ax.scatter(pt[~predriver_mask], adata.obs['supp_score'].values[~predriver_mask],
           c='#CCCCCC', s=1.5, alpha=0.3, linewidths=0)
ax.scatter(pt[predriver_mask], adata.obs['supp_score'].values[predriver_mask],
           c='#E64B35', s=3, alpha=0.8, linewidths=0, label='Pre-driver')
ax.set_xlabel('Pseudotime', fontsize=9)
ax.set_ylabel('Suppressive score', fontsize=9)
ax.set_title('Pre-driver candidates highlighted', fontsize=10, fontweight='bold')
ax.legend(fontsize=8)

fig.tight_layout()
for ext in ('png', 'pdf'):
    path = f'{FIG}/predriver_trajectory.{ext}'
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f'  Saved: {path}')
plt.close()

# ── Print summary ──────────────────────────────────────────────
print('\n' + '=' * 65)
print('Summary')
print('=' * 65)
print(f'  Myeloid cells analyzed : {adata.n_obs:,}')
print(f'  Datasets               : GSE125527 (UC) + GSE132465 (CRC)')
print(f'  Root cell subtype      : {adata.obs.loc[root_cell, "cell_subtype_fine"]}')
print(f'  Pseudotime range       : {pt.min():.3f} – {pt.max():.3f}')
print(f'  Terminals found        : {n_term}')
print(f'  Pre-driver candidates  : {n_pd:,} ({n_pd / adata.n_obs * 100:.1f}%)')
print(f'  Pre-driver pseudotime  : [{pt_lo:.3f}, {pt_hi:.3f}]')
print(f'\n  Results in : {RES}/')
print(f'  Figures in : {FIG}/')
print(f'\nDone at {datetime.now().strftime("%H:%M:%S")}')
