#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 2 — Fine sub-annotation v2 (hierarchical strategy)
T cells: CD8/CD4 lineage first → subtype (properly identifies Tex/Tpex)
Myeloid: adds TAM via FTH1/SOD2/SOCS3/OSM/HMOX1
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys, io, warnings
from datetime import datetime

warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

plt.rcParams['font.family']  = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

FIG_DIR = 'figures/phase2/fine_annotation'
RES_DIR = 'phase2_results/fine_annotation'
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

BATCH_UC = {'GSE125527', 'GSE150115'}


def savefig(fig, name):
    fig.savefig(f'{FIG_DIR}/{name}.png', dpi=200, bbox_inches='tight', facecolor='white')
    fig.savefig(f'{FIG_DIR}/{name}.pdf',           bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'   Saved: {FIG_DIR}/{name}.png/.pdf')


def normalize_copy(adata):
    a = adata.copy()
    a.X = a.layers['counts'].copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def score_module(adata_score, adata_target, genes, name):
    valid = [g for g in genes if g in adata_score.var_names]
    if len(valid) >= 2:
        sc.tl.score_genes(adata_score, valid, score_name=name, random_state=0)
        adata_target.obs[name] = adata_score.obs[name].values
    else:
        adata_target.obs[name] = 0.0


# ══════════════════════════════════════════════════════════════════
#  T CELL ANNOTATION — hierarchical
# ══════════════════════════════════════════════════════════════════

T_PALETTE = {
    'CD8_Teff':  '#E64B35',
    'CD8_Tex':   '#B22222',
    'CD8_Tpex':  '#FF8C69',
    'CD4_Tconv': '#4DBBD5',
    'Treg':      '#00A087',
    'Th17':      '#3C5488',
    'Tfh':       '#F39B7F',
    'gdT':       '#8491B4',
    'Cycling_T': '#91D1C2',
    'Unknown_T': '#CCCCCC',
}

T_DOTPLOT_GENES = [
    # CD8 lineage
    'CD8A', 'CD8B', 'GZMB', 'PRF1', 'NKG7', 'IFNG',       # Teff
    'HAVCR2', 'ENTPD1', 'TOX', 'LAG3', 'TIGIT', 'PDCD1',   # Tex
    'TCF7', 'SLAMF6', 'ID3', 'CXCR5',                       # Tpex
    'TRDC', 'TRGC1',                                          # gdT
    # CD4 lineage
    'CD4', 'IL7R', 'CCR7', 'SELL', 'KLF2',                  # Tconv
    'FOXP3', 'IL2RA', 'IKZF2',                               # Treg
    'CCR6', 'RORA', 'KLRB1',                                 # Th17
    'ICOS', 'BCL6', 'CD40LG',                                # Tfh
    # Cycling
    'STMN1', 'PCNA',
]


def annotate_tcells(adata_t):
    print('\n[T] Clustering (resolution=1.0) ...')
    sc.pp.neighbors(adata_t, use_rep='X_scvi', n_neighbors=30, random_state=0)
    sc.tl.leiden(adata_t, resolution=1.0, random_state=0, key_added='leiden')
    sc.tl.umap(adata_t, min_dist=0.3, random_state=0)
    print(f'    {adata_t.obs["leiden"].nunique()} clusters')

    print('[T] Scoring modules ...')
    a = normalize_copy(adata_t)

    # ── Lineage scores ─────────────────────────────────────────
    score_module(a, adata_t, ['CD8A', 'CD8B'],               'sc_CD8_lin')
    score_module(a, adata_t, ['CD4', 'IL7R'],                 'sc_CD4_lin')
    score_module(a, adata_t, ['TRDC', 'TRGC1', 'TRGC2'],     'sc_gdT_lin')
    score_module(a, adata_t, ['STMN1', 'PCNA'],               'sc_Cycling')

    # ── CD8 subtype scores ─────────────────────────────────────
    # Tex: exhaustion markers, penalise TCF7 via separate score
    score_module(a, adata_t, ['HAVCR2', 'ENTPD1', 'TOX', 'TIGIT', 'LAG3', 'PDCD1', 'CTLA4'],
                 'sc_Tex_pos')
    score_module(a, adata_t, ['TCF7', 'CCR7', 'SELL'],        'sc_Tex_neg')  # naive/progenitor → penalises Tex
    score_module(a, adata_t, ['TCF7', 'SLAMF6', 'ID3', 'CXCR5'],            'sc_Tpex')
    score_module(a, adata_t, ['GZMB', 'PRF1', 'IFNG', 'NKG7', 'GNLY', 'EOMES', 'TBX21'],
                 'sc_Teff')

    # Tex = pos score minus neg score
    adata_t.obs['sc_Tex'] = adata_t.obs['sc_Tex_pos'] - 0.5 * adata_t.obs['sc_Tex_neg']

    # ── CD4 subtype scores ─────────────────────────────────────
    score_module(a, adata_t, ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2'], 'sc_Treg')
    score_module(a, adata_t, ['CCR6', 'RORA', 'KLRB1'],             'sc_Th17')
    score_module(a, adata_t, ['CXCR5', 'ICOS', 'BCL6', 'CD40LG'],  'sc_Tfh')
    score_module(a, adata_t, ['TCF7', 'CCR7', 'SELL', 'KLF2', 'S1PR1', 'IL7R'],
                 'sc_Tconv')

    print('[T] Assigning subtypes per cluster (hierarchical) ...')
    cluster_map = {}
    for cl in sorted(adata_t.obs['leiden'].unique(), key=int):
        mask = adata_t.obs['leiden'] == cl
        s = adata_t.obs.loc[mask]
        n = mask.sum()

        # ① Cycling (highest priority)
        if s['sc_Cycling'].mean() > 0.6:
            cluster_map[cl] = 'Cycling_T'; continue

        cd8  = s['sc_CD8_lin'].mean()
        cd4  = s['sc_CD4_lin'].mean()
        gdt  = s['sc_gdT_lin'].mean()

        # ② gdT
        if gdt > 0.3 and gdt >= cd8 and gdt >= cd4:
            cluster_map[cl] = 'gdT'; continue

        # ③ CD8 lineage
        if cd8 >= cd4:
            scores = {
                'CD8_Tex':  s['sc_Tex'].mean(),
                'CD8_Tpex': s['sc_Tpex'].mean(),
                'CD8_Teff': s['sc_Teff'].mean(),
            }
            best = max(scores, key=scores.get)
            cluster_map[cl] = best if scores[best] > 0 else 'Unknown_T'

        # ④ CD4 lineage
        else:
            scores = {
                'Treg':      s['sc_Treg'].mean(),
                'Th17':      s['sc_Th17'].mean(),
                'Tfh':       s['sc_Tfh'].mean(),
                'CD4_Tconv': s['sc_Tconv'].mean(),
            }
            best = max(scores, key=scores.get)
            cluster_map[cl] = best if scores[best] > 0 else 'Unknown_T'

    adata_t.obs['subtype'] = adata_t.obs['leiden'].map(cluster_map).astype('category')

    print('\n    Cluster → subtype:')
    for cl, st in sorted(cluster_map.items(), key=lambda x: int(x[0])):
        n = (adata_t.obs['leiden'] == cl).sum()
        print(f'      {cl:>2s}: {st}  (n={n:,})')

    return adata_t


# ══════════════════════════════════════════════════════════════════
#  MYELOID ANNOTATION — adds TAM via FTH1/SOD2/SOCS3/OSM
# ══════════════════════════════════════════════════════════════════

MY_PALETTE = {
    'Classical_Mono':    '#E64B35',
    'NonClassical_Mono': '#F39B7F',
    'TAM':               '#B22222',
    'Inflam_Mac':        '#FF6600',
    'cDC1':              '#4DBBD5',
    'cDC2':              '#00A087',
    'pDC':               '#3C5488',
    'Proliferating_My':  '#91D1C2',
    'Unknown_My':        '#CCCCCC',
}

MY_DOTPLOT_GENES = [
    'CD14', 'S100A8', 'S100A9', 'FCN1', 'LYZ',          # Classical mono
    'FCGR3A', 'CX3CR1', 'CDKN1C',                         # NonClassical mono
    'CD68', 'FCGR1A', 'LGALS3',                            # Macrophage general
    'FTH1', 'SOD2', 'SOCS3', 'OSM', 'HMOX1',             # TAM (CRC fate)
    'IL1B', 'TNF',                                          # Inflam Mac
    'CADM1', 'IRF8', 'BATF3',                              # cDC1
    'CLEC10A', 'CD1C', 'FCER1A',                           # cDC2
    'IRF7', 'PLAC8', 'TCF4',                               # pDC
    'STMN1', 'PCNA',                                        # Proliferating
]


def annotate_myeloids(adata_my):
    print('\n[My] Clustering (resolution=0.8) ...')
    sc.pp.neighbors(adata_my, use_rep='X_scvi', n_neighbors=30, random_state=0)
    sc.tl.leiden(adata_my, resolution=0.8, random_state=0, key_added='leiden')
    sc.tl.umap(adata_my, min_dist=0.3, random_state=0)
    print(f'    {adata_my.obs["leiden"].nunique()} clusters')

    print('[My] Scoring modules ...')
    a = normalize_copy(adata_my)

    score_module(a, adata_my, ['CD14', 'S100A8', 'S100A9', 'FCN1', 'LYZ'],
                 'sc_ClassMono')
    score_module(a, adata_my, ['FCGR3A', 'CX3CR1', 'CDKN1C'],
                 'sc_NonClassMono')
    # TAM: CRC-fate macrophage markers (FTH1/SOD2 = oxidative stress/iron; OSM/SOCS3 = immunosuppression)
    score_module(a, adata_my, ['FTH1', 'SOD2', 'SOCS3', 'OSM', 'HMOX1', 'FCGR1A', 'LGALS3', 'THBS1', 'CD68'],
                 'sc_TAM')
    score_module(a, adata_my, ['IL1B', 'TNF', 'CD68', 'LYZ'],
                 'sc_InflamMac')
    score_module(a, adata_my, ['CADM1', 'IRF8', 'BATF3'],
                 'sc_cDC1')
    score_module(a, adata_my, ['CLEC10A', 'CD1C', 'FCER1A'],
                 'sc_cDC2')
    score_module(a, adata_my, ['IRF7', 'PLAC8', 'TCF4', 'GZMB'],
                 'sc_pDC')
    score_module(a, adata_my, ['STMN1', 'PCNA'],
                 'sc_Prolif')

    # Distinguish TAM from Inflam_Mac: TAM has low S100A8/FCN1 (non-monocyte-derived)
    score_module(a, adata_my, ['S100A8', 'S100A9', 'FCN1'],
                 'sc_mono_character')
    adata_my.obs['sc_TAM_adj'] = adata_my.obs['sc_TAM'] - 0.5 * adata_my.obs['sc_mono_character']

    print('[My] Assigning subtypes per cluster ...')
    cluster_map = {}
    for cl in sorted(adata_my.obs['leiden'].unique(), key=int):
        mask = adata_my.obs['leiden'] == cl
        s = adata_my.obs.loc[mask]

        if s['sc_Prolif'].mean() > 0.5:
            cluster_map[cl] = 'Proliferating_My'; continue

        scores = {
            'Classical_Mono':    s['sc_ClassMono'].mean(),
            'NonClassical_Mono': s['sc_NonClassMono'].mean(),
            'TAM':               s['sc_TAM_adj'].mean(),
            'Inflam_Mac':        s['sc_InflamMac'].mean(),
            'cDC1':              s['sc_cDC1'].mean(),
            'cDC2':              s['sc_cDC2'].mean(),
            'pDC':               s['sc_pDC'].mean(),
        }
        best = max(scores, key=scores.get)
        cluster_map[cl] = best if scores[best] > 0 else 'Unknown_My'

    adata_my.obs['subtype'] = adata_my.obs['leiden'].map(cluster_map).astype('category')

    print('\n    Cluster → subtype:')
    for cl, st in sorted(cluster_map.items(), key=lambda x: int(x[0])):
        n = (adata_my.obs['leiden'] == cl).sum()
        print(f'      {cl:>2s}: {st}  (n={n:,})')

    return adata_my


# ══════════════════════════════════════════════════════════════════
#  SHARED PLOTTING FUNCTIONS
# ══════════════════════════════════════════════════════════════════

def plot_umap_subtype(adata_sub, palette, title, name):
    umap = adata_sub.obsm['X_umap']
    shuf = np.random.default_rng(0).permutation(adata_sub.n_obs)
    S, A = 1.5, 0.5
    subtypes    = adata_sub.obs['subtype'].cat.categories.tolist()
    dis_colors  = {'UC': '#00A087', 'CRC': '#4DBBD5'}

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(title, fontsize=13, fontweight='bold')

    ax = axes[0]
    for st in subtypes:
        m = adata_sub.obs['subtype'].values == st
        ax.scatter(umap[m, 0], umap[m, 1], c=palette.get(st, '#999'),
                   s=S, alpha=A, linewidths=0, label=f'{st} ({m.sum():,})')
    ax.set_title('By subtype', fontsize=11, fontweight='bold')
    ax.axis('off')
    ax.legend(markerscale=5, fontsize=7, loc='upper left',
              bbox_to_anchor=(1.01, 1), borderaxespad=0)

    ax = axes[1]
    for dis, col in dis_colors.items():
        m = adata_sub.obs['disease'].values == dis
        ax.scatter(umap[shuf][m[shuf], 0], umap[shuf][m[shuf], 1],
                   c=col, s=S, alpha=A, linewidths=0, label=dis)
    ax.set_title('By disease', fontsize=11, fontweight='bold')
    ax.axis('off')
    ax.legend(markerscale=5, fontsize=9)

    plt.tight_layout()
    savefig(fig, name)


def plot_dotplot(adata_sub, gene_list, title, name):
    a = normalize_copy(adata_sub)
    genes = [g for g in gene_list if g in a.var_names]
    subtypes = adata_sub.obs['subtype'].cat.categories.tolist()

    means   = np.zeros((len(subtypes), len(genes)))
    pct_exp = np.zeros((len(subtypes), len(genes)))
    for i, st in enumerate(subtypes):
        mask = adata_sub.obs['subtype'].values == st
        X = a[mask, genes].X
        if hasattr(X, 'toarray'):
            X = X.toarray()
        means[i]   = X.mean(axis=0)
        pct_exp[i] = (X > 0).mean(axis=0) * 100

    means_norm = means / (means.max(axis=0, keepdims=True) + 1e-9)

    fig, ax = plt.subplots(figsize=(max(10, len(genes) * 0.38),
                                     max(4,  len(subtypes) * 0.5)))
    fig.suptitle(title, fontsize=11, fontweight='bold')

    for i, st in enumerate(subtypes):
        for j in range(len(genes)):
            size = (pct_exp[i, j] / 100) ** 0.5 * 220
            ax.scatter(j, i, s=size, c=[[plt.cm.Reds(means_norm[i, j])]],
                       linewidths=0.3, edgecolors='grey', alpha=0.9, zorder=2)

    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=7)
    ax.set_yticks(range(len(subtypes)))
    ax.set_yticklabels(subtypes, fontsize=9)
    ax.set_xlim(-0.5, len(genes) - 0.5)
    ax.set_ylim(-0.5, len(subtypes) - 0.5)
    ax.grid(axis='both', color='#EEEEEE', linewidth=0.5, zorder=1)
    ax.spines[['top', 'right']].set_visible(False)
    for pct, lbl in [(20, '20%'), (50, '50%'), (80, '80%')]:
        ax.scatter([], [], s=(pct/100)**0.5*220, c='grey', linewidths=0, label=lbl, alpha=0.7)
    ax.legend(title='% expressed', fontsize=7, title_fontsize=7,
              loc='upper left', bbox_to_anchor=(1.01, 1))
    plt.tight_layout()
    savefig(fig, name)


def plot_composition(adata_sub, palette, title, name):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(title, fontsize=11, fontweight='bold')

    for ax, col in zip(axes, ['disease', 'batch']):
        ct = pd.crosstab(adata_sub.obs[col], adata_sub.obs['subtype'],
                         normalize='index') * 100
        subtypes = ct.columns.tolist()
        bottom = np.zeros(len(ct))
        for st in subtypes:
            ax.bar(range(len(ct)), ct[st].values, bottom=bottom,
                   color=palette.get(st, '#999'), label=st, width=0.6)
            bottom += ct[st].values
        ax.set_xticks(range(len(ct)))
        ax.set_xticklabels(ct.index.tolist(), rotation=30, ha='right', fontsize=8)
        ax.set_ylabel('Proportion (%)', fontsize=9)
        ax.set_title(f'By {col}', fontsize=10, fontweight='bold')
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_ylim(0, 100)

    handles = [plt.Rectangle((0, 0), 1, 1, color=palette.get(st, '#999'))
               for st in ct.columns]
    axes[1].legend(handles, ct.columns.tolist(), fontsize=7,
                   loc='upper left', bbox_to_anchor=(1.01, 1))
    plt.tight_layout()
    savefig(fig, name)


def print_crosstab(adata_sub):
    print('\n    Subtype counts:')
    print(adata_sub.obs['subtype'].value_counts().to_string())
    print('\n    Subtype × disease:')
    print(pd.crosstab(adata_sub.obs['subtype'], adata_sub.obs['disease']).to_string())


# ══════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════
print('=' * 65)
print(f'Phase 2 — Fine annotation v2  |  {datetime.now().strftime("%H:%M:%S")}')
print('=' * 65)

print('\nLoading data ...')
adata_all = sc.read_h5ad('processed_data/merged_annotated.h5ad')
print(f'  Total: {adata_all.n_obs:,} cells')

if 'disease' not in adata_all.obs.columns:
    adata_all.obs['disease'] = adata_all.obs['batch'].map(
        lambda b: 'UC' if b in BATCH_UC else 'CRC')

adata_all.obs['cell_subtype_fine'] = adata_all.obs['cell_type_grouped'].astype(str)

# ── T cells ────────────────────────────────────────────────────
print('\n' + '=' * 65)
print(f'T cells  |  {datetime.now().strftime("%H:%M:%S")}')
print('=' * 65)
adata_t = adata_all[adata_all.obs['cell_type_grouped'] == 'T cells'].copy()
print(f'  {adata_t.n_obs:,} cells')

adata_t = annotate_tcells(adata_t)
print_crosstab(adata_t)

print('\n  Plotting ...')
plot_umap_subtype(adata_t, T_PALETTE, 'T cell fine annotation', 'tcell_umap_subtype')
plot_dotplot(adata_t, T_DOTPLOT_GENES, 'T cell marker genes', 'tcell_dotplot')
plot_composition(adata_t, T_PALETTE, 'T cell subtype composition', 'tcell_composition')

t_out = adata_t.obs[['batch', 'disease', 'cell_type_grouped', 'leiden', 'subtype']].copy()
t_out.columns = ['batch', 'disease', 'cell_type_grouped', 'leiden_T', 'subtype_T']
t_out.to_csv(f'{RES_DIR}/tcell_subtypes.csv')
print(f'  Saved: {RES_DIR}/tcell_subtypes.csv')

adata_all.obs.loc[adata_t.obs_names, 'cell_subtype_fine'] = adata_t.obs['subtype'].values

# ── Myeloid cells ───────────────────────────────────────────────
print('\n' + '=' * 65)
print(f'Myeloid cells  |  {datetime.now().strftime("%H:%M:%S")}')
print('=' * 65)
adata_my = adata_all[adata_all.obs['cell_type_grouped'] == 'Myeloids'].copy()
print(f'  {adata_my.n_obs:,} cells')

adata_my = annotate_myeloids(adata_my)
print_crosstab(adata_my)

print('\n  Plotting ...')
plot_umap_subtype(adata_my, MY_PALETTE, 'Myeloid cell fine annotation', 'myeloid_umap_subtype')
plot_dotplot(adata_my, MY_DOTPLOT_GENES, 'Myeloid marker genes', 'myeloid_dotplot')
plot_composition(adata_my, MY_PALETTE, 'Myeloid subtype composition', 'myeloid_composition')

my_out = adata_my.obs[['batch', 'disease', 'cell_type_grouped', 'leiden', 'subtype']].copy()
my_out.columns = ['batch', 'disease', 'cell_type_grouped', 'leiden_My', 'subtype_My']
my_out.to_csv(f'{RES_DIR}/myeloid_subtypes.csv')
print(f'  Saved: {RES_DIR}/myeloid_subtypes.csv')

adata_all.obs.loc[adata_my.obs_names, 'cell_subtype_fine'] = adata_my.obs['subtype'].values

# ── Save ───────────────────────────────────────────────────────
print('\n' + '=' * 65)
print('Saving merged_annotated_fine.h5ad ...')
adata_all.write_h5ad('processed_data/merged_annotated_fine.h5ad')

print('\nFinal cell_subtype_fine:')
print(adata_all.obs['cell_subtype_fine'].value_counts().to_string())
print('\n' + '=' * 65)
print(f'Done  |  {datetime.now().strftime("%H:%M:%S")}')
print('=' * 65)
