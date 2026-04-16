#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2 – Panel D: InfiltrationScore barplot
Compares Pre-driver / SPP1+ fraction at Normal / Border / Tumor
per patient and as group mean, demonstrating the border-enrichment
"physical barrier" hypothesis.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import os

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
plt.rcParams['font.family'] = 'Arial'

RES = 'phase2_results/spatial_v2'
FIG_OUT = 'figures/phase2/spatial_v2'
os.makedirs(FIG_OUT, exist_ok=True)

# ── Load infiltration scores ──────────────────────────────────────
df = pd.read_csv(f'{RES}/infiltration_scores_v2.csv')
print(df.to_string())

CLASS_ORDER  = ['Normal', 'Border', 'Tumor']
CLASS_COLORS = {'Normal': '#4DBBD5', 'Border': '#E64B35', 'Tumor': '#F39B7F'}

metrics = [
    ('frac_SPP1',         'SPP1+ fraction',         'SPP1+ Myeloid Fraction'),
    ('frac_suppressive',  'Suppressive fraction',    'Suppressive Myeloid Fraction'),
    ('mean_Predriver_v2', 'Pre-driver score (mean)', 'Pre-driver Score'),
]

fig, axes = plt.subplots(1, 3, figsize=(12, 5))
fig.subplots_adjust(left=0.08, right=0.97, top=0.88, bottom=0.12, wspace=0.38)

for ax, (col, ylabel, title) in zip(axes, metrics):
    means, sems, patients = {}, {}, {}
    for cls in CLASS_ORDER:
        vals = df.loc[df['Class'] == cls, col].values
        means[cls]    = vals.mean()
        sems[cls]     = vals.std() / np.sqrt(len(vals))
        patients[cls] = df.loc[df['Class'] == cls, ['Patient', col]]

    x = np.arange(len(CLASS_ORDER))
    bars = ax.bar(x, [means[c] for c in CLASS_ORDER],
                  color=[CLASS_COLORS[c] for c in CLASS_ORDER],
                  width=0.55, alpha=0.85,
                  yerr=[sems[c] for c in CLASS_ORDER],
                  capsize=4, error_kw={'linewidth': 1.2})

    # Jitter patient dots
    rng = np.random.default_rng(42)
    for xi, cls in enumerate(CLASS_ORDER):
        pts = patients[cls][col].values
        jx  = xi + rng.uniform(-0.18, 0.18, size=len(pts))
        ax.scatter(jx, pts, c='black', s=28, alpha=0.7,
                   linewidths=0.5, edgecolors='white', zorder=5)

    # Significance brackets: Normal vs Border, Normal vs Tumor
    ymax = max(df[col].max() * 1.05, 0.01)
    pairs = [('Normal', 'Border'), ('Normal', 'Tumor')]
    bracket_y = ymax * 1.08
    for (g1, g2) in pairs:
        v1 = df.loc[df['Class'] == g1, col].values
        v2 = df.loc[df['Class'] == g2, col].values
        _, pval = stats.mannwhitneyu(v1, v2, alternative='two-sided')
        stars = '***' if pval < 0.001 else ('**' if pval < 0.01 else ('*' if pval < 0.05 else 'ns'))
        x1, x2 = CLASS_ORDER.index(g1), CLASS_ORDER.index(g2)
        h = bracket_y
        ax.plot([x1, x1, x2, x2], [h, h + ymax * 0.03, h + ymax * 0.03, h],
                lw=1.0, c='#444')
        ax.text((x1 + x2) / 2, h + ymax * 0.04, stars, ha='center', fontsize=9)
        bracket_y += ymax * 0.14

    ax.set_xticks(x)
    ax.set_xticklabels(CLASS_ORDER, fontsize=9.5)
    ax.set_ylabel(ylabel, fontsize=9.5)
    ax.set_title(title, fontsize=10, pad=4)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=8.5)

    if col == 'frac_SPP1':
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda y, _: f'{y*100:.0f}%'))

axes[0].text(-0.20, 1.04, 'D', transform=axes[0].transAxes,
             fontsize=14, fontweight='bold', va='bottom')

fig.suptitle('Infiltration Score — Myeloid Suppressive Enrichment by Region (GSE144735)',
             fontsize=10.5, y=0.97)

for ext in ('png', 'pdf'):
    path = f'{FIG_OUT}/infiltration_score_panelD.{ext}'
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f'Saved: {path}')
plt.close()
print('Done.')
