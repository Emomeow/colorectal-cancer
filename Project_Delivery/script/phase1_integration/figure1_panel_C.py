#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 1 Panel C — Cell type composition per dataset (stacked bar)
"""
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os, sys, io, shutil

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

os.makedirs('figures/phase1', exist_ok=True)

print("Loading data ...")
adata = sc.read_h5ad('processed_data/merged_annotated.h5ad')
print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

# ── Colour palette — same as Panel A/B/D (Nature npg) ────────────────────────
TYPE_COLORS = {
    'T cells'   : '#E64B35',
    'NK cells'  : '#3C5488',
    'B cells'   : '#4DBBD5',
    'Myeloids'  : '#00A087',
    'Mast cells': '#F39B7F',
    'Epithelial': '#91D1C2',
    'Stromal'   : '#DC0000',
    'Malignant' : '#B09C85',
    'Unknown'   : '#AAAAAA',
}

CT_ORDER = [
    'T cells', 'NK cells', 'B cells', 'Myeloids', 'Mast cells',
    'Epithelial', 'Stromal', 'Malignant', 'Unknown',
]

# ── Batch order: CRC first, then UC ──────────────────────────────────────────
BATCH_ORDER = ['GSE132257', 'GSE132465', 'cellxgene', 'GSE125527', 'GSE150115']
batch_labels = []
for b in BATCH_ORDER:
    n = (adata.obs['batch'] == b).sum()
    batch_labels.append(f"{b}\n({n:,} cells)")

# ── Compute proportions ───────────────────────────────────────────────────────
ct = pd.crosstab(adata.obs['batch'], adata.obs['cell_type_grouped'])
for col in CT_ORDER:
    if col not in ct.columns:
        ct[col] = 0
ct   = ct[CT_ORDER].loc[BATCH_ORDER]
prop = ct.div(ct.sum(axis=1), axis=0)

print("\nComposition (proportions):")
print(prop.round(3).to_string())

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 5.5))

bottoms = np.zeros(len(BATCH_ORDER))
x = np.arange(len(BATCH_ORDER))

for ct_name in CT_ORDER:
    vals = prop[ct_name].values
    ax.bar(x, vals, bottom=bottoms,
           color=TYPE_COLORS[ct_name], width=0.65,
           edgecolor='white', linewidth=0.5, label=ct_name)
    bottoms += vals

# ── Axis formatting ───────────────────────────────────────────────────────────
ax.set_xticks(x)
ax.set_xticklabels(batch_labels, fontsize=9)
ax.set_ylim(0, 1)
ax.set_ylabel('Proportion of cells', fontsize=10)
ax.set_xlabel('Dataset', fontsize=10, labelpad=8)
ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1, decimals=0))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='y', labelsize=9)

# ── Disease label brackets — placed BELOW x-tick labels ──────────────────────
# Use blended transform: x = data coords, y = axes fraction (negative = below)
trans = ax.get_xaxis_transform()

for (lo, hi, label) in [(0, 2, 'CRC'), (3, 4, 'UC')]:
    xlo = lo - 0.38
    xhi = hi + 0.38
    xmid = (xlo + xhi) / 2
    y_line = -0.22   # just below x-tick labels (axes fraction)
    y_text = -0.285  # below the bracket line

    # Horizontal bracket line
    ax.annotate('', xy=(xhi, y_line), xytext=(xlo, y_line),
                arrowprops=dict(arrowstyle='-', color='#666666', lw=1.0),
                xycoords=('data', 'axes fraction'),
                textcoords=('data', 'axes fraction'),
                annotation_clip=False)
    # Short vertical ticks at ends
    for xend in [xlo, xhi]:
        ax.annotate('', xy=(xend, y_line - 0.015),
                    xytext=(xend, y_line),
                    arrowprops=dict(arrowstyle='-', color='#666666', lw=1.0),
                    xycoords=('data', 'axes fraction'),
                    textcoords=('data', 'axes fraction'),
                    annotation_clip=False)
    # Label text below bracket
    ax.text(xmid, y_text, label, ha='center', va='top',
            fontsize=9, fontweight='bold', color='#333333',
            transform=trans, clip_on=False)

# ── Title — placed above axes, safely above any bracket ──────────────────────
ax.set_title('Figure 1C — Cell Type Composition per Dataset',
             fontsize=13, fontweight='bold', pad=14)

# ── Legend ────────────────────────────────────────────────────────────────────
handles = [mpatches.Patch(facecolor=TYPE_COLORS[c], label=c, edgecolor='#cccccc', linewidth=0.5)
           for c in CT_ORDER]
ax.legend(handles=handles, title='Cell type', title_fontsize=9,
          fontsize=8, bbox_to_anchor=(1.01, 1), loc='upper left',
          framealpha=0.85, edgecolor='#cccccc', handlelength=1.2)

plt.tight_layout()

# ── Save ──────────────────────────────────────────────────────────────────────
out_png = 'figures/phase1/figure1_panel_C_composition.png'
out_pdf = 'figures/phase1/figure1_panel_C_composition.pdf'
plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(out_pdf, bbox_inches='tight', facecolor='white')
plt.close()
print(f"\nSaved: {out_png}")
print(f"Saved: {out_pdf}")
print("Done.")
