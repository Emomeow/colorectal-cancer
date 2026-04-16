#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd, numpy as np, sys, io, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

plt.rcParams['font.family']  = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

TBL = 'Project_Delivery/results/Figure2_Pre_driver_Trajectory/Spatial_Distribution/tables/'
OUT = 'Project_Delivery/results/Figure2_Pre_driver_Trajectory/Spatial_Distribution/figures/'

df_sub  = pd.read_csv(TBL + 'border_enrichment_subtypes.csv')
df_gene = pd.read_csv(TBL + 'border_enrichment_curated_genes.csv')
df_top  = pd.read_csv(TBL + 'border_enrichment_top200_genes_fullscan.csv')

CT_COLOR = {
    'Myeloids':        '#D62728',
    'T cells':         '#1F77B4',
    'B cells':         '#FF7F0E',
    'Stromal cells':   '#2CA02C',
    'Mast cells':      '#9467BD',
    'Epithelial cells':'#CCCCCC',
}

fig, axes = plt.subplots(1, 3, figsize=(16, 7),
                          gridspec_kw={'width_ratios': [1.2, 1, 1.4]})
fig.suptitle('Tumor Border Enrichment Analysis (GSE144735, 6 patients)',
             fontsize=12, fontweight='bold', y=1.01)

# ── Panel A: Cell subtype lollipop ────────────────────────────────
ax = axes[0]
ax.set_title('A  Cell Subtypes', fontsize=10, fontweight='bold', loc='left')

sub_plot = df_sub[df_sub['log2FC_Border_Normal'].notna()].copy()
sub_plot = sub_plot[sub_plot['Cell_subtype'] != 'Epithelial cells']
sub_plot = sub_plot.sort_values('log2FC_Border_Normal').reset_index(drop=True)

for i, row in sub_plot.iterrows():
    col   = CT_COLOR.get(row['Cell_type'], '#888888')
    fc    = row['log2FC_Border_Normal']
    fdr   = row['FDR_Border_vs_Normal']
    alpha = 1.0 if fdr < 0.1 else 0.4
    lw    = 2.0 if fdr < 0.1 else 0.8
    ax.plot([0, fc], [i, i], color=col, lw=lw, alpha=alpha, solid_capstyle='round')
    ms = 10 if fdr < 0.05 else (7 if fdr < 0.1 else 5)
    ax.scatter(fc, i, s=ms**2, color=col, alpha=alpha,
               edgecolors='black' if fdr < 0.05 else col, linewidths=0.8, zorder=5)
    if fdr < 0.05:
        ax.text(fc + 0.12, i, '**', fontsize=8, va='center', color='#333333')
    elif fdr < 0.1:
        ax.text(fc + 0.12, i, '*', fontsize=8, va='center', color='#333333')

ax.set_yticks(range(len(sub_plot)))
ax.set_yticklabels(sub_plot['Cell_subtype'].tolist(), fontsize=8)
ax.axvline(0, color='#888888', lw=0.8, ls='--')
ax.set_xlabel('log2FC (Border / Normal)', fontsize=9)
ax.spines[['top', 'right']].set_visible(False)
ax.grid(axis='x', color='#EEEEEE', lw=0.5)

legend_els = [mpatches.Patch(facecolor=v, label=k)
              for k, v in CT_COLOR.items() if k != 'Epithelial cells']
ax.legend(handles=legend_els, fontsize=7, loc='lower right',
          frameon=True, framealpha=0.8)

# ── Panel B: Curated gene bubble plot ─────────────────────────────
ax = axes[1]
ax.set_title('B  Curated Immune Genes', fontsize=10, fontweight='bold', loc='left')

gene_plot = df_gene[df_gene['log2FC_Border_Normal'] > 0].copy()
gene_plot = gene_plot.sort_values('log2FC_Border_Normal').reset_index(drop=True)

for i, row in gene_plot.iterrows():
    fc    = row['log2FC_Border_Normal']
    fdr   = row['FDR']
    col   = CT_COLOR.get(row['Main_cell_type'], '#888888')
    alpha = 1.0 if fdr < 0.1 else 0.4
    lw    = 2.0 if fdr < 0.1 else 0.8
    size  = max(20, -np.log10(fdr + 1e-3) * 40)
    ax.plot([0, fc], [i, i], color=col, lw=lw, alpha=alpha, solid_capstyle='round')
    ax.scatter(fc, i, s=size, color=col, alpha=alpha,
               edgecolors='black' if fdr < 0.1 else col, linewidths=0.7, zorder=5)
    if fdr < 0.1:
        label = '**' if fdr < 0.05 else '*'
        ax.text(fc + 0.06, i, label, fontsize=8, va='center', color='#333333')

ax.set_yticks(range(len(gene_plot)))
ax.set_yticklabels(gene_plot['Gene'].tolist(), fontsize=8.5, fontstyle='italic')
ax.axvline(0, color='#888888', lw=0.8, ls='--')
ax.set_xlabel('log2FC (Border / Normal)', fontsize=9)
ax.spines[['top', 'right']].set_visible(False)
ax.grid(axis='x', color='#EEEEEE', lw=0.5)

for fdr_v, lbl in [(0.05, 'FDR=0.05'), (0.1, 'FDR=0.1'), (0.3, 'FDR=0.3')]:
    s = max(20, -np.log10(fdr_v + 1e-3) * 40)
    ax.scatter([], [], s=s, color='grey', alpha=0.6, label=lbl)
ax.legend(title='-log10(FDR)', fontsize=7, title_fontsize=7,
          loc='lower right', frameon=True, framealpha=0.8)

# ── Panel C: Top genes from full scan grouped by cell type ────────
ax = axes[2]
ax.set_title('C  Top Border Genes (Full Scan, log2FC>2)', fontsize=10,
             fontweight='bold', loc='left')

top_plot = df_top[df_top['log2FC'] > 2.0].copy()
group_order = ['Myeloids', 'Stromal cells', 'T cells', 'B cells']
top_plot['ct_rank'] = top_plot['Main_cell_type'].map(
    {ct: i for i, ct in enumerate(group_order)}).fillna(99)
top_plot = (top_plot.sort_values(['ct_rank', 'log2FC'], ascending=[True, False])
                    .groupby('Main_cell_type').head(8)
                    .sort_values(['ct_rank', 'log2FC'], ascending=[True, False])
                    .reset_index(drop=True))

y_positions, y_labels, y_colors, y_fdrs, y_fcs = [], [], [], [], []
y_offset = 0
prev_ct  = None
separator_ys = []
group_ys = {}

for _, row in top_plot.iterrows():
    ct = row['Main_cell_type']
    if ct != prev_ct and prev_ct is not None:
        separator_ys.append(y_offset - 0.3)
        y_offset += 0.6
    if ct not in group_ys:
        group_ys[ct] = []
    group_ys[ct].append(y_offset)
    y_positions.append(y_offset)
    y_labels.append(row['Gene'])
    y_colors.append(CT_COLOR.get(ct, '#888888'))
    y_fdrs.append(row['FDR'])
    y_fcs.append(row['log2FC'])
    prev_ct = ct
    y_offset += 1

for sep_y in separator_ys:
    ax.axhline(sep_y, color='#DDDDDD', lw=1.0)

for y, fc, col, fdr in zip(y_positions, y_fcs, y_colors, y_fdrs):
    alpha = 1.0 if fdr < 0.15 else 0.45
    lw    = 1.8 if fdr < 0.15 else 0.8
    ax.plot([0, fc], [y, y], color=col, lw=lw, alpha=alpha, solid_capstyle='round')
    ax.scatter(fc, y, s=36, color=col, alpha=alpha,
               edgecolors='black' if fdr < 0.15 else col,
               linewidths=0.6, zorder=5)

ax.set_yticks(y_positions)
ax.set_yticklabels(y_labels, fontsize=7.5, fontstyle='italic')
ax.axvline(0, color='#888888', lw=0.8, ls='--')
ax.set_xlabel('log2FC (Border / Normal)', fontsize=9)
ax.spines[['top', 'right']].set_visible(False)
ax.grid(axis='x', color='#EEEEEE', lw=0.5)

x_max = max(y_fcs) * 1.02
for ct, ys in group_ys.items():
    mid_y = (min(ys) + max(ys)) / 2
    col   = CT_COLOR.get(ct, '#888888')
    ax.text(x_max, mid_y, ct, fontsize=7, color=col,
            fontweight='bold', va='center', ha='left')

ax.set_xlim(right=x_max + 1.5)

# ── Shared footnote ───────────────────────────────────────────────
fig.text(0.5, -0.01,
         '** FDR<0.05   * FDR<0.1   . FDR<0.2   ns not significant  |  '
         'Mann-Whitney U (per-patient fractions), BH correction',
         ha='center', fontsize=7.5, color='#555555')

plt.tight_layout()
fname = 'border_enrichment_summary'
fig.savefig(f'{OUT}{fname}.png', dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(f'{OUT}{fname}.pdf',           bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {OUT}{fname}.png / .pdf')
