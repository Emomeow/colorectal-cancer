#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Novel TF Candidates Figure — Phase 2
Two-panel figure:
  Panel A: CRC TAM vs UC Myeloid pseudobulk DEG — novel TF lollipop
           (MAF/BHLHE41/BATF3 up; CIITA/NR4A1/IRF8/SREBF2/ZEB2 down)
  Panel B: PPARG & CEBPB decoupleR delta activity (TAM vs Mono)
           + Pseudotime correlation for existing SCENIC TFs
"""
import pandas as pd, numpy as np, sys, io, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

plt.rcParams['font.family']  = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

OUT = 'Project_Delivery/results/Figure2_Pre_driver_Trajectory/SCENIC_Regulon/figures/'
os.makedirs(OUT, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────
deg = pd.read_csv('phase2_results/pseudobulk/deg_TAM.csv')
delta = pd.read_csv('phase2_results/regulon/delta_tf_TAM_vs_Mono.csv',
                    header=0, names=['TF', 'delta_activity'])
pt_corr = pd.read_csv('phase2_results/tf_v2/tf_pseudotime_correlation.csv')

# Novel TF selection: UP and DOWN in CRC TAM (vs UC Myeloids, FDR<0.05)
# Focused on non-textbook, mechanistically interesting
TF_UP   = ['MAF', 'BATF3', 'BHLHE41', 'HES4']
TF_DOWN = ['CIITA', 'NR4A1', 'IRF8', 'NR4A2', 'CEBPD', 'SREBF2', 'ZEB2', 'TBX21']
TF_ORDER = TF_UP[::-1] + TF_DOWN  # down then up for lollipop (bottom=down)

TF_ANNOT = {
    'MAF':    'c-MAF: IL-10/TGFB1 suppressive program',
    'BATF3':  'DC-like suppressive phenotype',
    'BHLHE41':'Tissue-resident macrophage identity (DEC2)',
    'HES4':   'Notch target; inflammatory attenuation',
    'CIITA':  'MHC-II master regulator: antigen presentation LOSS',
    'NR4A1':  'Nur77: homeostatic Mc; loss drives pro-tumor TAM',
    'IRF8':   'DC1/M1 TF; type-I IFN capacity loss',
    'NR4A2':  'NR4A family; co-loss with NR4A1',
    'CEBPD':  'CEBP delta; inflammatory differentiation TF loss',
    'SREBF2': 'De novo lipid synthesis TF (dn); exogenous uptake (OLR1) (up)',
    'ZEB2':   'EMT/myeloid differentiation regulator',
    'TBX21':  'T-bet: M1 pro-inflammatory; lost in CRC TAM',
}

t_sel = deg[deg['gene'].isin(TF_ORDER)].set_index('gene').reindex(TF_ORDER)
t_sel = t_sel.dropna(subset=['log2FC'])

fig, axes = plt.subplots(1, 2, figsize=(16, 7),
                          gridspec_kw={'width_ratios': [1.8, 1.0]})
fig.suptitle('Novel Transcriptional Regulators in CRC TAM\n'
             '(CRC TAM vs UC Myeloids, Pseudobulk DEG)',
             fontsize=11, fontweight='bold', y=1.01)

# ── Panel A: Novel TF lollipop ────────────────────────────────────────
ax = axes[0]
ax.set_title('A  CRC-specific TF Gain / Loss in TAM', fontsize=10,
             fontweight='bold', loc='left')

valid = [g for g in TF_ORDER if g in t_sel.index]
n = len(valid)

UP_COLOR   = '#E64B35'   # red for up
DOWN_COLOR = '#3C5488'   # blue for down
STAR_TFS   = {'MAF', 'CIITA', 'NR4A1', 'PPARG'}  # top-story TFs

for i, gene in enumerate(valid):
    fc  = t_sel.loc[gene, 'log2FC']
    fdr = t_sel.loc[gene, 'fdr']
    col = UP_COLOR if fc > 0 else DOWN_COLOR
    is_star = gene in STAR_TFS

    alpha = 1.0 if fdr < 0.001 else (0.8 if fdr < 0.01 else 0.6)
    lw    = 2.5 if is_star else 1.8
    ms    = 100 if is_star else (64 if fdr < 0.001 else 36)

    ax.plot([0, fc], [i, i], color=col, lw=lw, alpha=alpha,
            solid_capstyle='round')
    ax.scatter(fc, i, s=ms, color=col, alpha=alpha,
               edgecolors='black', linewidths=0.8, zorder=5)

    # FDR label
    if fdr < 1e-20:  sig = '***'
    elif fdr < 1e-10: sig = '**'
    elif fdr < 0.001: sig = '*'
    else:             sig = ''
    x_off = 0.12 if fc >= 0 else -0.12
    ha    = 'left' if fc >= 0 else 'right'
    if sig:
        ax.text(fc + x_off, i, sig, fontsize=7.5, va='center',
                color='#333333', ha=ha)

    # Star marker for top-story TFs
    if is_star:
        side = fc + 0.28 if fc >= 0 else fc - 0.28
        ha2  = 'left' if fc >= 0 else 'right'
        ax.text(side, i, '[*]', fontsize=9, va='center',
                color='#B8860B', ha=ha2)

ax.set_yticks(range(n))
ax.set_yticklabels(valid, fontsize=9, fontstyle='italic')
for lbl, gene in zip(ax.get_yticklabels(), valid):
    fc = t_sel.loc[gene, 'log2FC']
    col = UP_COLOR if fc > 0 else DOWN_COLOR
    if gene in STAR_TFS:
        lbl.set_fontweight('bold')
    lbl.set_color(col)

# Separator between up and down
sep_y = len(TF_DOWN) - 0.5
ax.axhline(sep_y, color='#CCCCCC', lw=1.2, ls='--')
ax.text(-0.05, sep_y + 0.3, 'UP in CRC TAM (up)', fontsize=7.5,
        color=UP_COLOR, ha='right', va='bottom', fontstyle='italic')
ax.text(-0.05, sep_y - 0.3, 'DOWN in CRC TAM (dn)', fontsize=7.5,
        color=DOWN_COLOR, ha='right', va='top', fontstyle='italic')

ax.axvline(0, color='#888888', lw=0.8, ls='--')
ax.set_xlabel('log2FC (CRC TAM / UC Myeloids)', fontsize=9)
ax.spines[['top', 'right']].set_visible(False)
ax.grid(axis='x', color='#EEEEEE', lw=0.5)

# Annotation box on right side
for i, gene in enumerate(valid):
    if gene in TF_ANNOT:
        fc = t_sel.loc[gene, 'log2FC']
        col = UP_COLOR if fc > 0 else DOWN_COLOR
        # Only annotate key genes
        if gene in ('MAF', 'CIITA', 'NR4A1', 'BHLHE41', 'SREBF2', 'IRF8'):
            ax.annotate(TF_ANNOT[gene],
                        xy=(fc, i),
                        xytext=(fc + (0.8 if fc > 0 else -0.8), i),
                        fontsize=6.5, color='#444444', va='center',
                        ha='left' if fc > 0 else 'right',
                        arrowprops=dict(arrowstyle='-', color='#AAAAAA', lw=0.5))

# ── Panel B: Delta activity (decoupleR) + Pseudotime correlation ──────
ax = axes[1]
ax.set_title('B  TF Activity\n(decoupleR: TAM vs Mono)', fontsize=10,
             fontweight='bold', loc='left')

# Plot all delta activity TFs
DELTA_COLOR = '#E64B35'   # positive delta = UP in TAM
delta_sorted = delta.sort_values('delta_activity')

# Highlight story TFs
STORY_DELTA = {'PPARG', 'CEBPB', 'IRF1', 'RELA', 'NFKB1', 'MYC',
               'NR4A1', 'BATF'}

for i, row in delta_sorted.reset_index(drop=True).iterrows():
    tf  = row['TF']
    val = row['delta_activity']
    col = '#E64B35' if val > 0 else '#3C5488'
    is_story = tf in STORY_DELTA
    is_lipid = tf in ('PPARG', 'CEBPB')

    alpha = 1.0 if is_story else 0.55
    lw    = 2.2 if is_lipid else (1.8 if is_story else 1.0)
    ms    = 100 if is_lipid else (64 if is_story else 25)

    ax.plot([0, val], [i, i], color=col, lw=lw, alpha=alpha,
            solid_capstyle='round')
    ax.scatter(val, i, s=ms, color=col, alpha=alpha,
               edgecolors='black' if is_story else col,
               linewidths=0.7, zorder=5)

    if is_lipid:
        ax.text(val + (0.15 if val >= 0 else -0.15), i, '[*]',
                fontsize=9, va='center',
                color='#B8860B', ha='left' if val >= 0 else 'right')

ax.set_yticks(range(len(delta_sorted)))
ax.set_yticklabels(delta_sorted['TF'].tolist(), fontsize=8)
for lbl, (_, row) in zip(ax.get_yticklabels(),
                          delta_sorted.reset_index(drop=True).iterrows()):
    tf  = row['TF']
    val = row['delta_activity']
    col = '#E64B35' if val > 0 else '#3C5488'
    if tf in STORY_DELTA:
        lbl.set_fontweight('bold')
    if tf in ('PPARG', 'CEBPB'):
        lbl.set_color('#2CA02C')
    else:
        lbl.set_color(col)

ax.axvline(0, color='#888888', lw=0.8, ls='--')
ax.set_xlabel('D Activity (TAM - Classical_Mono)', fontsize=9)
ax.spines[['top', 'right']].set_visible(False)
ax.grid(axis='x', color='#EEEEEE', lw=0.5)

# Annotations for PPARG and CEBPB
for _, row in delta_sorted.iterrows():
    if row['TF'] in ('PPARG', 'CEBPB'):
        i = delta_sorted.reset_index(drop=True).index[
            delta_sorted.reset_index(drop=True)['TF'] == row['TF']].tolist()
        if i:
            label = 'Lipid TF\n(PPARG)' if row['TF'] == 'PPARG' else 'M2 TF\n(CEBPB)'
            ax.annotate(label,
                        xy=(row['delta_activity'], i[0]),
                        xytext=(row['delta_activity'] + 0.5, i[0]),
                        fontsize=6.5, color='#2CA02C', va='center',
                        arrowprops=dict(arrowstyle='->', color='#2CA02C', lw=0.8))

# ── Legend ────────────────────────────────────────────────────────────
legend_els = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=UP_COLOR,
           markersize=8, label='UP in CRC TAM'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=DOWN_COLOR,
           markersize=8, label='DOWN in CRC TAM'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#2CA02C',
           markersize=8, label='Lipid metabolism TF'),
    Line2D([0], [0], marker='*', color='w', markerfacecolor='#B8860B',
           markersize=10, label='Story TF ([*])'),
]
axes[0].legend(handles=legend_els, fontsize=7.5, loc='lower right',
               frameon=True, framealpha=0.85)

fig.text(0.5, -0.01,
         'Panel A: CRC TAM vs UC Myeloid pseudobulk DEG (scVI batch-corrected, Wilcoxon, BH-FDR)  |  '
         '*** FDR<10-\u00b2\u2070  ** FDR<10-\u00b9\u2070  * FDR<0.001  '
         '[*] mechanistically novel, CRC-specific\n'
         'Panel B: decoupleR-ULM TF activity (CollecTRI network); DActivity = TAM mean - Classical_Mono mean  |  '
         'Green labels: lipid metabolism TFs supporting OLR1+ TAM model',
         ha='center', fontsize=7, color='#555555')

plt.tight_layout()
fname = 'tf_novel_candidates'
fig.savefig(f'{OUT}{fname}.png', dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(f'{OUT}{fname}.pdf',           bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {OUT}{fname}.png / .pdf')
