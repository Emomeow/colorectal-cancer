#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Border Enrichment Summary v2
Highlights OLR1+ lipid-associated TAM + FAP+ Myofibroblast dual-barrier story
Panel A: Cell subtypes (SPP1+B renamed OLR1+ TAM, Myofibroblast starred)
Panel B: OLR1+ TAM lipid signature + Myofibroblast matrix genes (story-driven)
Panel C: Full-scan top genes by cell type (OLR1/MMP9/TREM2 highlighted)
"""
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

# Novel story genes: OLR1+ TAM lipid signature + Myofibroblast markers
# Pull from full-scan table (more complete log2FC data)
STORY_GENES_MY  = ['OLR1', 'TREM2', 'TREM1', 'MMP9', 'SPP1', 'SCD', 'APOC1']
STORY_GENES_STR = ['FAP', 'INHBA', 'PODNL1', 'COL10A1', 'MMP1', 'SUGCT', 'UNC5B']

def get_gene_row(gene, df_top, df_curated):
    """Fetch log2FC/FDR for a gene, prefer full-scan, fallback to curated."""
    row = df_top[df_top['Gene'] == gene]
    if not row.empty:
        r = row.iloc[0]
        return {'Gene': gene,
                'Main_cell_type': r['Main_cell_type'],
                'log2FC_Border_Normal': r['log2FC'],
                'FDR': r['FDR']}
    row = df_curated[df_curated['Gene'] == gene]
    if not row.empty:
        r = row.iloc[0]
        return {'Gene': gene,
                'Main_cell_type': r['Main_cell_type'],
                'log2FC_Border_Normal': r.get('log2FC_Border_Normal', r.get('log2FC', 0)),
                'FDR': r['FDR']}
    return None

story_rows = []
for g in STORY_GENES_MY + STORY_GENES_STR:
    r = get_gene_row(g, df_top, df_gene)
    if r:
        story_rows.append(r)
df_story = pd.DataFrame(story_rows)

fig, axes = plt.subplots(1, 3, figsize=(17, 7),
                          gridspec_kw={'width_ratios': [1.3, 1.1, 1.4]})
fig.suptitle('Tumor Border Enrichment Analysis\n'
             'OLR1+ TAM + FAP+ Myofibroblast Dual-Barrier Axis (GSE144735, 6 patients)',
             fontsize=11, fontweight='bold', y=1.01)

# ── Panel A: Cell subtype lollipop ────────────────────────────────────
ax = axes[0]
ax.set_title('A  Cell Subtypes', fontsize=10, fontweight='bold', loc='left')

sub_plot = df_sub[df_sub['log2FC_Border_Normal'].notna()].copy()
sub_plot = sub_plot[sub_plot['Cell_subtype'] != 'Epithelial cells']
sub_plot = sub_plot.sort_values('log2FC_Border_Normal').reset_index(drop=True)

# Rename SPP1+B → OLR1+ TAM (SPP1+B) in display
RENAME_MAP = {'SPP1+B': 'OLR1+ TAM (SPP1+B)', 'SPP1+A': 'SPP1+A TAM'}
HIGHLIGHT_CELLS = {'OLR1+ TAM (SPP1+B)', 'Myofibroblasts'}

for i, row in sub_plot.iterrows():
    col   = CT_COLOR.get(row['Cell_type'], '#888888')
    fc    = row['log2FC_Border_Normal']
    fdr   = row['FDR_Border_vs_Normal']
    alpha = 1.0 if fdr < 0.1 else 0.4
    lw    = 2.0 if fdr < 0.1 else 0.8
    disp_name = RENAME_MAP.get(row['Cell_subtype'], row['Cell_subtype'])
    is_highlight = disp_name in HIGHLIGHT_CELLS or row['Cell_subtype'] == 'Myofibroblasts'

    if is_highlight:
        ax.plot([0, fc], [i, i], color=col, lw=2.5, alpha=1.0,
                solid_capstyle='round')
        ax.scatter(fc, i, s=144, color=col, alpha=1.0,
                   edgecolors='black', linewidths=1.2, zorder=5)
    else:
        ax.plot([0, fc], [i, i], color=col, lw=lw, alpha=alpha,
                solid_capstyle='round')
        ms = 10 if fdr < 0.05 else (7 if fdr < 0.1 else 5)
        ax.scatter(fc, i, s=ms**2, color=col, alpha=alpha,
                   edgecolors='black' if fdr < 0.05 else col,
                   linewidths=0.8, zorder=5)

    if fdr < 0.05:
        ax.text(fc + 0.15, i, '**', fontsize=8, va='center', color='#333333')
    elif fdr < 0.1:
        ax.text(fc + 0.15, i, '*', fontsize=8, va='center', color='#333333')

    if is_highlight:
        star = ' [*]' if row['Cell_subtype'] in ('SPP1+B', 'Myofibroblasts') else ''
        ax.text(fc + 0.15, i, star, fontsize=9, va='center', color='#B8860B')

display_labels = [RENAME_MAP.get(r, r) for r in sub_plot['Cell_subtype'].tolist()]
ax.set_yticks(range(len(sub_plot)))
ax.set_yticklabels(display_labels, fontsize=7.5)
# Bold highlight labels
for lbl in ax.get_yticklabels():
    if lbl.get_text() in HIGHLIGHT_CELLS or 'OLR1' in lbl.get_text() or 'Myofibroblast' in lbl.get_text():
        lbl.set_fontweight('bold')
        lbl.set_color('#333333')

ax.axvline(0, color='#888888', lw=0.8, ls='--')
ax.set_xlabel('log2FC (Border / Normal)', fontsize=9)
ax.spines[['top', 'right']].set_visible(False)
ax.grid(axis='x', color='#EEEEEE', lw=0.5)

legend_els = [mpatches.Patch(facecolor=v, label=k)
              for k, v in CT_COLOR.items() if k != 'Epithelial cells']
ax.legend(handles=legend_els, fontsize=7, loc='lower right',
          frameon=True, framealpha=0.8)

# ── Panel B: Story genes lollipop ─────────────────────────────────────
ax = axes[1]
ax.set_title('B  OLR1+ TAM + Myofibroblast\nSignature Genes', fontsize=9,
             fontweight='bold', loc='left')

# Separate myeloid vs stromal; myeloid first (sorted), then separator, then stromal
df_my_s  = df_story[df_story['Main_cell_type'] == 'Myeloids'].sort_values(
    'log2FC_Border_Normal').reset_index(drop=True)
df_str_s = df_story[df_story['Main_cell_type'] == 'Stromal cells'].sort_values(
    'log2FC_Border_Normal').reset_index(drop=True)

y_offset = 0
y_positions, y_labels, y_fcs, y_fdrs, y_cols = [], [], [], [], []

for _, row in df_my_s.iterrows():
    y_positions.append(y_offset); y_labels.append(row['Gene'])
    y_fcs.append(row['log2FC_Border_Normal']); y_fdrs.append(row['FDR'])
    y_cols.append(CT_COLOR['Myeloids']); y_offset += 1

sep_y = y_offset - 0.3
y_offset += 0.8

for _, row in df_str_s.iterrows():
    y_positions.append(y_offset); y_labels.append(row['Gene'])
    y_fcs.append(row['log2FC_Border_Normal']); y_fdrs.append(row['FDR'])
    y_cols.append(CT_COLOR['Stromal cells']); y_offset += 1

ax.axhline(sep_y, color='#CCCCCC', lw=1.0, ls=':')

OLR1_HIGHLIGHT = {'OLR1', 'TREM2', 'TREM1', 'MMP9'}
FAP_HIGHLIGHT  = {'FAP', 'INHBA', 'COL10A1', 'MMP1'}

for y, fc, col, fdr, gene in zip(y_positions, y_fcs, y_cols, y_fdrs, y_labels):
    is_key = gene in OLR1_HIGHLIGHT or gene in FAP_HIGHLIGHT
    alpha = 1.0 if (fdr < 0.2 or is_key) else 0.45
    lw    = 2.0 if is_key else (1.5 if fdr < 0.2 else 0.8)
    size  = max(20, -np.log10(fdr + 1e-3) * 35)
    ax.plot([0, fc], [y, y], color=col, lw=lw, alpha=alpha,
            solid_capstyle='round')
    ax.scatter(fc, y, s=size, color=col, alpha=alpha,
               edgecolors='black' if fdr < 0.15 else col,
               linewidths=0.7, zorder=5)
    if fdr < 0.15:
        lbl = '**' if fdr < 0.05 else '*'
        ax.text(fc + 0.12, y, lbl, fontsize=8, va='center', color='#333333')

ax.set_yticks(y_positions)
ax.set_yticklabels(y_labels, fontsize=8.5, fontstyle='italic')
for lbl in ax.get_yticklabels():
    g = lbl.get_text()
    if g in OLR1_HIGHLIGHT:
        lbl.set_color(CT_COLOR['Myeloids']); lbl.set_fontweight('bold')
    elif g in FAP_HIGHLIGHT:
        lbl.set_color(CT_COLOR['Stromal cells']); lbl.set_fontweight('bold')

ax.axvline(0, color='#888888', lw=0.8, ls='--')
ax.set_xlabel('log2FC (Border / Normal)', fontsize=9)
ax.spines[['top', 'right']].set_visible(False)
ax.grid(axis='x', color='#EEEEEE', lw=0.5)

# Group labels on right
x_max_b = max(y_fcs) * 1.02 if max(y_fcs) > 0 else 2
for group, color, rows in [
    ('OLR1+ TAM', CT_COLOR['Myeloids'], df_my_s),
    ('Myofibroblast', CT_COLOR['Stromal cells'], df_str_s)
]:
    idxs = list(range(len(df_my_s))) if 'OLR1' in group else \
           list(range(len(df_my_s)+1, len(y_positions)))
    if idxs:
        mid = (y_positions[idxs[0]] + y_positions[idxs[-1]]) / 2
        ax.text(x_max_b + 0.1, mid, group, fontsize=7, color=color,
                fontweight='bold', va='center', ha='left')
ax.set_xlim(right=x_max_b + 1.8)

for fdr_v, lbl in [(0.05, 'FDR=0.05'), (0.1, 'FDR=0.1'), (0.2, 'FDR=0.2')]:
    s = max(20, -np.log10(fdr_v + 1e-3) * 35)
    ax.scatter([], [], s=s, color='grey', alpha=0.6, label=lbl)
ax.legend(title='-log10(FDR)', fontsize=7, title_fontsize=7,
          loc='lower right', frameon=True, framealpha=0.8)

# ── Panel C: Top genes grouped by cell type (full scan) ───────────────
ax = axes[2]
ax.set_title('C  Top Border Genes (Full Scan, log2FC>1.5)', fontsize=10,
             fontweight='bold', loc='left')

top_plot = df_top[df_top['log2FC'] > 1.5].copy()
group_order = ['Myeloids', 'Stromal cells', 'T cells', 'B cells']
top_plot['ct_rank'] = top_plot['Main_cell_type'].map(
    {ct: i for i, ct in enumerate(group_order)}).fillna(99)
top_plot = (top_plot.sort_values(['ct_rank', 'log2FC'], ascending=[True, False])
                    .groupby('Main_cell_type').head(8)
                    .sort_values(['ct_rank', 'log2FC'], ascending=[True, False])
                    .reset_index(drop=True))

PANEL_C_HIGHLIGHT = {'OLR1', 'MMP9', 'TREM2', 'TREM1', 'FAP', 'INHBA', 'COL10A1', 'MMP1'}

y_positions_c, y_labels_c, y_colors_c = [], [], []
y_fdrs_c, y_fcs_c, y_highlight_c = [], [], []
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
    y_positions_c.append(y_offset)
    y_labels_c.append(row['Gene'])
    y_colors_c.append(CT_COLOR.get(ct, '#888888'))
    y_fdrs_c.append(row['FDR'])
    y_fcs_c.append(row['log2FC'])
    y_highlight_c.append(row['Gene'] in PANEL_C_HIGHLIGHT)
    prev_ct = ct
    y_offset += 1

for sep_y in separator_ys:
    ax.axhline(sep_y, color='#DDDDDD', lw=1.0)

for y, fc, col, fdr, highlight in zip(y_positions_c, y_fcs_c, y_colors_c,
                                       y_fdrs_c, y_highlight_c):
    alpha = 1.0 if fdr < 0.15 else 0.45
    lw    = 2.2 if highlight else (1.8 if fdr < 0.15 else 0.8)
    ax.plot([0, fc], [y, y], color=col, lw=lw, alpha=alpha,
            solid_capstyle='round')
    ms = 64 if highlight else 36
    ax.scatter(fc, y, s=ms, color=col, alpha=alpha,
               edgecolors='black' if (fdr < 0.15 or highlight) else col,
               linewidths=0.6, zorder=5)

ax.set_yticks(y_positions_c)
ax.set_yticklabels(y_labels_c, fontsize=7.5, fontstyle='italic')
for lbl, hl, col in zip(ax.get_yticklabels(), y_highlight_c, y_colors_c):
    if hl:
        lbl.set_fontweight('bold')
        lbl.set_color(col)

ax.axvline(0, color='#888888', lw=0.8, ls='--')
ax.set_xlabel('log2FC (Border / Normal)', fontsize=9)
ax.spines[['top', 'right']].set_visible(False)
ax.grid(axis='x', color='#EEEEEE', lw=0.5)

x_max = max(y_fcs_c) * 1.02
for ct, ys in group_ys.items():
    mid_y = (min(ys) + max(ys)) / 2
    col   = CT_COLOR.get(ct, '#888888')
    ax.text(x_max, mid_y, ct, fontsize=7, color=col,
            fontweight='bold', va='center', ha='left')
ax.set_xlim(right=x_max + 1.6)

# ── Footnote ──────────────────────────────────────────────────────────
fig.text(0.5, -0.01,
         '** FDR<0.05   * FDR<0.1   [*] novel story gene  |  '
         'Mann-Whitney U (per-patient fractions), BH correction  |  '
         'Bold: OLR1+ TAM lipid signature (red) or FAP+ Myofibroblast matrix genes (green)',
         ha='center', fontsize=7, color='#555555')

plt.tight_layout()
fname = 'border_enrichment_summary_v2'
fig.savefig(f'{OUT}{fname}.png', dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(f'{OUT}{fname}.pdf',           bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {OUT}{fname}.png / .pdf')
